import argparse
import datetime
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional

import multiprocessing as mp
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--trace_dir', help='path to trace directory', required=True)
parser.add_argument('--results_dir', help='path to results directory', required=True)
parser.add_argument(
    '--submissions_dir',
    help='path to Submissions directory (default: ./Submissions)',
    default='./Submissions',
)
parser.add_argument(
    '--jobs',
    type=int,
    default=0,
    help='number of parallel trace runs (0 = use all cores)',
)
parser.add_argument(
    '--rebuild_binaries',
    action='store_true',
    help='force rebuild of per-predictor cbp binaries even if they already exist',
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _resolve_repo_relative(path_str: str) -> Path:
    p = Path(path_str).expanduser()
    if p.is_absolute():
        return p
    return (REPO_ROOT / p).resolve()


# Map trace-group/workload -> which submission predictor to use.
# Based on reports/predictor_best_by_trace_group.txt (winner by mean CycWPPKI).
WORKLOAD_TO_SUBMISSION = {
    'compress': 'REGISTER-VALUE-AWARE-TORU-KOIZUMI',
    'fp': 'REGISTER-VALUE-AWARE-TORU-KOIZUMI',
    'infra': 'LOAD-VALUE-CORRELATOR-MAN',
    'int': 'REGISTER-VALUE-AWARE-TORU-KOIZUMI',
    'media': 'MULTI-PERSPECTIVE-PREDICTOR-JIMENEZ',
    'web': 'REGISTER-VALUE-AWARE-TORU-KOIZUMI',
}


def _safe_name(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('_')


def choose_submission_for_trace_path(trace_path: Path) -> str:
    # Expected trace layout: .../<workload>/<run>_trace.gz
    workload = trace_path.parent.name
    return WORKLOAD_TO_SUBMISSION.get(workload, WORKLOAD_TO_SUBMISSION['int'])


def build_cbps_for_submissions(
    *,
    needed_submissions: set[str],
    bins_root: Path,
    submissions_dir: Path,
    force_rebuild: bool,
) -> dict[str, Path]:
    """Build a cbp binary for each submission in an isolated build dir.

    Returns: submission_name -> cbp_binary_path
    """
    if not submissions_dir.exists():
        raise FileNotFoundError(f"Submissions directory not found: {submissions_dir}")

    bins_root.mkdir(parents=True, exist_ok=True)
    build_root = bins_root / '_oracle_build'
    build_root.mkdir(parents=True, exist_ok=True)

    # We compile only the simulator shim objects; lib is shared from the repo.
    core_files = [
        'Makefile',
        'cbp.h',
        'cond_branch_predictor_interface.cc',
        'cond_branch_predictor_interface.h',
        'my_cond_branch_predictor.cc',
        'my_cond_branch_predictor.h',
        'cbp2016_tage_sc_l.h',
        'cbp2016_tage_sc_l_192kb.h',
    ]

    out: dict[str, Path] = {}
    for submission_name in sorted(needed_submissions):
        sub_dir = submissions_dir / submission_name
        if not sub_dir.is_dir():
            raise FileNotFoundError(f"Submission not found: {sub_dir}")

        bin_dir = bins_root / _safe_name(submission_name)
        cbp_out = bin_dir / 'cbp'
        if cbp_out.exists() and not force_rebuild:
            out[submission_name] = cbp_out
            continue

        # Fresh isolated build dir
        build_dir = build_root / _safe_name(submission_name)
        if build_dir.exists():
            shutil.rmtree(build_dir)
        build_dir.mkdir(parents=True, exist_ok=True)

        # Provide lib/ (prefer the submission's lib/ if present; otherwise use repo lib/).
        submission_lib = sub_dir / 'lib'
        if submission_lib.is_dir():
            shutil.copytree(submission_lib, build_dir / 'lib')
        else:
            # Symlink to avoid copying large/static artifacts.
            os.symlink(REPO_ROOT / 'lib', build_dir / 'lib')

        # Copy baseline core files (submission will overwrite as needed).
        for rel in core_files:
            src = REPO_ROOT / rel
            if src.exists():
                shutil.copy2(src, build_dir / rel)

        # Overlay all submission files into the build dir.
        for child in sub_dir.iterdir():
            if child.name == '.DS_Store':
                continue
            if child.name == 'lib':
                continue
            dst = build_dir / child.name
            if child.is_dir():
                shutil.copytree(child, dst)
            else:
                shutil.copy2(child, dst)

        # Build
        jobs = mp.cpu_count()
        subprocess.run(
            ['make', 'clean'],
            cwd=str(build_dir),
            check=True,
        )
        subprocess.run(
            ['make', f'-j{jobs}'],
            cwd=str(build_dir),
            check=True,
        )

        if not (build_dir / 'cbp').exists():
            raise RuntimeError(f"Build succeeded but cbp missing: {(build_dir / 'cbp')}")

        bin_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(build_dir / 'cbp', cbp_out)
        # Ensure executable bit is preserved/available.
        cbp_out.chmod(0o755)
        out[submission_name] = cbp_out

    return out


# Worker globals (set via Pool initializer; important on macOS where multiprocessing uses spawn)
WORKER_RESULTS_DIR: Optional[Path] = None
WORKER_SUBMISSION_TO_CBP: Optional[dict[str, Path]] = None


def _init_worker(results_dir_str: str, submission_to_cbp: dict[str, Path]) -> None:
    global WORKER_RESULTS_DIR, WORKER_SUBMISSION_TO_CBP
    WORKER_RESULTS_DIR = Path(results_dir_str)
    WORKER_SUBMISSION_TO_CBP = submission_to_cbp

def get_trace_paths(start_path):
    ret_list = []
    for root, dirs, files in os.walk(start_path):
        for my_file in files:
            if(my_file.endswith('_trace.gz')):
                ret_list.append(os.path.join(root, my_file))
    return ret_list

def process_run_op(pass_status, my_trace_path, my_run_name, op_file):
    run_name_split = re.split(r"\/", my_run_name)
    wl_name = run_name_split[0]
    run_name = run_name_split[1]
    print(f'Extracting data from : {op_file} |  WL:{wl_name} | Run:{run_name}')
    exec_time = 0

    _Instr = 0
    _Cycles = 0
    _IPC = 0
    _NumBr = 0
    _MispBr = 0
    _BrPerCyc = 0
    _MispBrPerCyc = 0
    _MR = 0
    _MPKI = 0
    _CycWP = 0
    _CycWPAvg = 0
    _CycWPPKI = 0
    
    _50PercInstr = 0
    _50PercCycles = 0
    _50PercIPC = 0
    _50PercNumBr = 0
    _50PercMispBr = 0
    _50PercBrPerCyc = 0
    _50PercMispBrPerCyc = 0
    _50PercMR = 0
    _50PercMPKI = 0
    _50PercCycWP = 0
    _50PercCycWPAvg = 0
    _50PercCycWPPKI = 0

    trace_size = os.path.getsize(my_trace_path) / (1024 * 1024)

    pass_status_str = 'Fail'

    process_50perc_section = False
    process_100perc_section = False

    _50perc_section_header = 'DIRECT CONDITIONAL BRANCH PREDICTION MEASUREMENTS (50 Perc instructions)'
    _100perc_section_header = 'DIRECT CONDITIONAL BRANCH PREDICTION MEASUREMENTS (Full Simulation i.e. Counts Not Reset When Warmup Ends)'

    found_50perc_line_to_process = False
    found_100perc_line_to_process = False

    if pass_status:
        pass_status_str = 'Pass'
        with open(op_file, "r") as text_file:
            #for line in my_run_output.splitlines():
            for line in text_file:
                if not line.strip():
                    continue

                if('ExecTime' in line):
                    exec_time = line.strip().split()[-1]

                if(not process_50perc_section and _50perc_section_header in line):
                    process_50perc_section = True
                    process_100perc_section = False

                if(not process_100perc_section and _100perc_section_header in line):
                    process_50perc_section = False
                    process_100perc_section = True


                if process_50perc_section:
                    if (found_50perc_line_to_process):
                        curr_split_line = line.split()
                        _50PercInstr              = curr_split_line[0]
                        _50PercCycles             = curr_split_line[1]
                        _50PercIPC                = curr_split_line[2]
                        _50PercNumBr              = curr_split_line[3]
                        _50PercMispBr             = curr_split_line[4]
                        _50PercBrPerCyc           = curr_split_line[5]
                        _50PercMispBrPerCyc       = curr_split_line[6]
                        _50PercMR                 = curr_split_line[7]
                        _50PercMPKI               = curr_split_line[8]
                        _50PercCycWP              = curr_split_line[9]
                        _50PercCycWPAvg           = curr_split_line[10]
                        _50PercCycWPPKI           = curr_split_line[11]
                        process_50perc_section = False
                        found_50perc_line_to_process = False

                    if(all(x in line for x in ['Instr', 'Cycles', 'IPC', 'NumBr', 'MispBr', 'BrPerCyc', 'MispBrPerCyc', 'MR', 'MPKI', 'CycWP', 'CycWPAvg'])):
                        found_50perc_line_to_process = True

                if process_100perc_section:
                    if (found_100perc_line_to_process):
                        curr_split_line = line.split()
                        _Instr              = curr_split_line[0]
                        _Cycles             = curr_split_line[1]
                        _IPC                = curr_split_line[2]
                        _NumBr              = curr_split_line[3]
                        _MispBr             = curr_split_line[4]
                        _BrPerCyc           = curr_split_line[5]
                        _MispBrPerCyc       = curr_split_line[6]
                        _MR                 = curr_split_line[7]
                        _MPKI               = curr_split_line[8]
                        _CycWP              = curr_split_line[9]
                        _CycWPAvg           = curr_split_line[10]
                        _CycWPPKI           = curr_split_line[11]
                        process_100perc_section = False
                        found_100perc_line_to_process = False
                    if(all(x in line for x in ['Instr', 'Cycles', 'IPC', 'NumBr', 'MispBr', 'BrPerCyc', 'MispBrPerCyc', 'MR', 'MPKI', 'CycWP', 'CycWPAvg'])):
                        found_100perc_line_to_process = True

    retval = {
            'Workload'                : wl_name,
            'Run'                     : run_name,
            'TraceSize'               : trace_size,
            'Status'                  : pass_status_str,
            'ExecTime'                : exec_time,

            'Instr'                   : _Instr,
            'Cycles'                  : _Cycles,
            'IPC'                     : _IPC,
            'NumBr'                   : _NumBr,
            'MispBr'                  : _MispBr,
            'BrPerCyc'                : _BrPerCyc,
            'MispBrPerCyc'            : _MispBrPerCyc,
            'MR'                      : _MR,
            'MPKI'                    : _MPKI,
            'CycWP'                   : _CycWP,
            'CycWPAvg'                : _CycWPAvg,
            'CycWPPKI'                : _CycWPPKI,

            '50PercInstr'             : _50PercInstr,
            '50PercCycles'            : _50PercCycles,
            '50PercIPC'               : _50PercIPC,
            '50PercNumBr'             : _50PercNumBr,
            '50PercMispBr'            : _50PercMispBr,
            '50PercBrPerCyc'          : _50PercBrPerCyc,
            '50PercMispBrPerCyc'      : _50PercMispBrPerCyc,
            '50PercMR'                : _50PercMR,
            '50PercMPKI'              : _50PercMPKI,
            '50PercCycWP'             : _50PercCycWP,
            '50PercCycWPAvg'          : _50PercCycWPAvg,
            '50PercCycWPPKI'          : _50PercCycWPPKI,
    }
    return retval

def execute_trace(my_trace_path):
    assert(os.path.exists(my_trace_path))

    if WORKER_RESULTS_DIR is None or WORKER_SUBMISSION_TO_CBP is None:
        raise RuntimeError('Worker not initialized (missing results dir or cbp binaries)')
        
    run_split = re.split(r"\/", my_trace_path)
    my_wl = run_split[-2] 
    # traces/int/int_0_trace.gz
    run_name = run_split[-1].split(".")[-2]
    if not os.path.exists(f'{WORKER_RESULTS_DIR}/{my_wl}'):
        if not os.path.exists(f'{WORKER_RESULTS_DIR}/{my_wl}'):
            os.makedirs(f'{WORKER_RESULTS_DIR}/{my_wl}', exist_ok=True)

    do_process = True
    my_run_name = f'{my_wl}/{run_name}'
    submission_name = choose_submission_for_trace_path(Path(my_trace_path))
    cbp_bin = WORKER_SUBMISSION_TO_CBP[submission_name]
    exec_cmd = [str(cbp_bin), str(my_trace_path)]
    op_file = f'{WORKER_RESULTS_DIR}/{my_wl}/{run_name}__{_safe_name(submission_name)}.log'
    if os.path.exists(op_file):
        #print(f"OP file:{op_file} already exists. Not running again!")
        do_process = False
    #
    pass_status = True
    if do_process:
        print(f'Begin processing run:{my_run_name}')
        try:
            begin_time = time.time()
            run_op = subprocess.check_output(exec_cmd, text=True)
            end_time = time.time()
            exec_time = end_time - begin_time
            with open(op_file, "w") as text_file:
                print(f"CMD:{' '.join(exec_cmd)}", file=text_file)
                print(f"PredictorSubmission:{submission_name}", file=text_file)
                print(f"{run_op}", file=text_file)
                print(f"ExecTime = {exec_time}", file=text_file)
        except:
            print(f'Run: {my_run_name} failed')
            pass_status = False
    return(pass_status, my_trace_path, op_file, my_run_name)



if __name__ == '__main__':
    args = parser.parse_args()
    trace_dir = _resolve_repo_relative(args.trace_dir)
    results_dir = _resolve_repo_relative(args.results_dir)
    submissions_dir = _resolve_repo_relative(args.submissions_dir)

    my_traces = get_trace_paths(trace_dir)
    print(f'Got {len(my_traces)} traces')

    results_dir.mkdir(parents=True, exist_ok=True)

    if len(my_traces) == 0:
        df = pd.DataFrame(
            columns=[
                'Workload',
                'Run',
                'TraceSize',
                'Status',
                'ExecTime',
                'Instr',
                'Cycles',
                'IPC',
                'NumBr',
                'MispBr',
                'BrPerCyc',
                'MispBrPerCyc',
                'MR',
                'MPKI',
                'CycWP',
                'CycWPAvg',
                'CycWPPKI',
                '50PercInstr',
                '50PercCycles',
                '50PercIPC',
                '50PercNumBr',
                '50PercMispBr',
                '50PercBrPerCyc',
                '50PercMispBrPerCyc',
                '50PercMR',
                '50PercMPKI',
                '50PercCycWP',
                '50PercCycWPAvg',
                '50PercCycWPPKI',
            ]
        )
        df.to_csv(results_dir / 'results.csv', index=False)
        print(f"No traces found under: {trace_dir}")
        print(f"Wrote empty results.csv to: {results_dir / 'results.csv'}")
        raise SystemExit(0)

    # Pre-build the predictor binaries needed for this trace set.
    needed_submissions = {choose_submission_for_trace_path(Path(p)) for p in my_traces}
    bins_root = results_dir / '_oracle_bins'
    submission_to_cbp = build_cbps_for_submissions(
        needed_submissions=needed_submissions,
        bins_root=bins_root,
        submissions_dir=submissions_dir,
        force_rebuild=bool(args.rebuild_binaries),
    )

    # For parallel runs using all CPU cores:
    num_cores = args.jobs if args.jobs and args.jobs > 0 else mp.cpu_count()
    total_traces = len(my_traces)
    print(f'\n=== Using {num_cores} CPU cores for parallel trace execution ===')
    print(f'=== Total traces to process: {total_traces} ===\n')

    results = []
    completed = 0
    with mp.Pool(
        processes=num_cores,
        initializer=_init_worker,
        initargs=(str(results_dir), submission_to_cbp),
    ) as pool:
        # Use imap_unordered for progress tracking
        for result in pool.imap_unordered(execute_trace, my_traces):
            results.append(result)
            completed += 1
            print(f'Progress: {completed}/{total_traces} traces completed ({100*completed//total_traces}%)', flush=True)

    print(f'\n=== All {total_traces} traces completed! ===\n')

    # For serial runs:
    #results = []
    #for my_trace in my_traces:
    #    results.append(execute_trace(my_trace))
    
    df = pd.DataFrame(columns=['Workload', 'Run', 'TraceSize', 'ExecTime', 'Instr', 'Cycles', 'IPC', 'NumBr', 'MispBr', 'BrPerCyc', 'MispBrPerCyc', 'MR', 'MPKI', 'CycWP',  'CycWPAvg', 'CycWPPKI', '50PercInstr', '50PercCycles', '50PercIPC', '50PercNumBr', '50PercMispBr', '50PercBrPerCyc', '50PercMispBrPerCyc', '50PercMR', '50PercMPKI', '50PercCycWP', '50PercCycWPAvg', '50PercCycWPPKI'])
    for my_result in results:
        pass_status = my_result[0]
        trace_path = my_result[1]
        op_file = my_result[2]
        my_run_name = my_result[3]
        run_dict = {}
        run_dict = process_run_op(pass_status, trace_path, my_run_name, op_file)
        my_df = pd.DataFrame([run_dict])
        if not df.empty:
            df = pd.concat([df, my_df], ignore_index=True)
        else:
            df = my_df.copy()
    print(df)
    df.to_csv(f'{results_dir}/results.csv', index=False)
    
    
    unique_wls = df['Workload'].unique()
    
    print('\n\n----------------------------------Aggregate Metrics Per Workload Category----------------------------------\n')
    for my_wl in unique_wls:
        my_wl_br_misp_pki_amean = df[df['Workload'] == my_wl]['50PercMPKI'].astype(float).mean()
        my_wl_cyc_wp_pki_amean = df[df['Workload'] == my_wl]['50PercCycWPPKI'].astype(float).mean()
        print(f'WL:{my_wl:<10} Branch Misprediction PKI(BrMisPKI) AMean : {my_wl_br_misp_pki_amean}')
        print(f'WL:{my_wl:<10} Cycles On Wrong-Path PKI(CycWpPKI) AMean : {my_wl_cyc_wp_pki_amean}')
    print('-----------------------------------------------------------------------------------------------------------')
    
    br_misp_pki_amean = df['50PercMPKI'].astype(float).mean()
    cyc_wp_pki_amean = df['50PercCycWPPKI'].astype(float).mean()
    #ipc_geomean = df['50PercIPC'].astype(float).apply(gmean)
    
    print('\n\n---------------------------------------------Aggregate Metrics---------------------------------------------\n')
    print(f'Branch Misprediction PKI(BrMisPKI) AMean : {br_misp_pki_amean}')
    print(f'Cycles On Wrong-Path PKI(CycWpPKI) AMean : {cyc_wp_pki_amean}')
    print('-----------------------------------------------------------------------------------------------------------')
