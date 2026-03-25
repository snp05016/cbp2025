# CBP2025 Graph + Report Generation 
Here’s how to get graphs and reports out of this repo without having to remember where every script lives.

also: there are *two* different “pipelines”, and they output to different places.

1) **Single-predictor analysis graphs** (one predictor at a time) using `generateReportForPredictor.sh`.
2) **Multi-predictor comparison graphs** (comparing many predictors at once) using scripts in `scripts/analysis/`.

---

## 1) Single-predictor graphs (one predictor, one `results.csv`)

### Where `results.csv` needs to be

The scripts that `runPythonScripts.sh` calls live in:

- `scripts/analysis/report_generators/`

They all read a file named **`results.csv`** from the **current working directory**.
so if you `cd` somewhere else and run them directly, make sure `results.csv` is in that directory.

The driver `runPythonScripts.sh` tries to be helpful: it searches for `results.csv` in a few common places:

- repo root: `./results.csv`
- `results/results.csv`
- `scripts/analysis/report_generators/results.csv`

Then it copies whatever it found into `scripts/analysis/report_generators/results.csv` (if needed) before running the Python scripts.

### Required columns (minimum-ish)

Most plots assume these columns exist in `results.csv`:

- Identity / grouping:
  - `Run` (e.g. `int_11_trace`)
  - `Workload` (category label used by some plots)
- Accuracy + difficulty:
  - `MPKI`, `50PercMPKI`
  - `MR` (usually a percent string like `"12.34%"`)
  - `50PercMR` (usually percent string)
- Performance:
  - `IPC`, `50PercIPC`
  - `Cycles`, `50PercCycles`
  - `Instr`
- Branch density / counts:
  - `BrPerCyc`, `NumBr`, `MispBr`
- Wrong-path cost:
  - `CycWPAvg`, `CycWPPKI`

If a plot crashes, it’s usually because one of these columns is missing or spelled differently, or because a “number” is stored as a string.

### What generates graphs (and where they end up)

All of the scripts below are invoked by `runPythonScripts.sh` (it does a `cd scripts/analysis/report_generators` first).
graphs now always land under the **repo-root** `reports/` directory.

- `scripts/analysis/report_generators/generate_graphs.py`
  - Writes many `.png` files to: `reports/<branch_predictor_name>/04_performance_metrics/graphs/`
  - Examples: `category_avg_misprediction_rate.png`, `*_performance_metrics.png`, `mpki_boxplot_benchmarks.png`, …

- `scripts/analysis/report_generators/graphComparator.py`
  - Writes: `reports/<branch_predictor_name>/01_misprediction_analysis/graphs/<branch_predictor_name>__mpki_comparison_Representative_Mix.png`

- `scripts/analysis/report_generators/compareMPKI_vs50PercMPKI.py`
  - Writes: `reports/<branch_predictor_name>/03_phase_behavior/graphs/<branch_predictor_name>__mpki_consistency_<category>.png`
  - Note: the script defaults to `category = 'int'` inside the file; change it if you want `fp`, `web`, etc.

- `scripts/analysis/report_generators/compareMetrics.py`
  - Writes: `reports/<branch_predictor_name>/05_comparative_analysis/divergent_behavior/<branch_predictor_name>__discrepancy_<category>_<different>_vs_<similar>.png`
  - Note: defaults inside the file (e.g. `category='fp'`, metrics) control what it emits.

- `scripts/analysis/report_generators/rangeOfBranchesAcrossBenchmarks.py`
  - Writes: `reports/<branch_predictor_name>/02_branch_density/graphs/<branch_predictor_name>__branch_density_geomean.png`

Text-only analysis (no graphs) also run by the pipeline:

- `scripts/analysis/report_generators/baselineDifficulty.py`
  - Writes text to the file that `runPythonScripts.sh` redirects to: `reports/<branch_predictor_name>/06_difficulty_analysis/intrinsic_hardness/baseline_difficulty_report.txt`

### How to run the single-predictor pipeline

From repo root (predictor name is required). you can copy/paste these:

```bash
chmod +x generateReportForPredictor.sh
./generateReportForPredictor.sh <branch_predictor_name> [optional/path/to/csv]

# examples
./generateReportForPredictor.sh baseline
./generateReportForPredictor.sh tage-sc-l-alberto-ros results/tage-sc-l-alberto-ros-results.csv
```

`runPythonScripts.sh` still works too (it’s just a wrapper that calls `generateReportForPredictor.sh`).

Outputs:

- Root report tree: `reports/<branch_predictor_name>/`
- Graphs are saved inside that folder and also prefixed with `<branch_predictor_name>__`.

Example:
- `reports/baseline/04_performance_metrics/graphs/baseline__mpki_boxplot_benchmarks.png`

---

## 2) Multi-predictor comparison graphs (many CSVs at once)

These scripts are not part of `runPythonScripts.sh` by default, but they generate the “compare predictors against each other” plots and tables.

### Where the predictor CSVs should live

This repo already has a `results/` folder containing multiple predictor result CSVs, e.g.:

- `results/baseline.csv`
- `results/programming-idiom-predictor-mose.csv`
- `results/tage-sc-l-alberto-ros-results.csv`
- `results/tage-scl-andrez-seznec-results.csv`
- `results/load-value-correlator-man-results.csv`
- `results/register-value-aware-toru-results.csv`

### The scripts you actually run

If you want a single command that runs both steps (per-predictor report *and* multi-predictor comparisons), use:

```bash
chmod +x runPredictorAndCompare.sh
./runPredictorAndCompare.sh <branch_predictor_name> [optional/path/to/csv]
```

- `scripts/analysis/compare_predictors.py`
  - Automatically discovers all `*.csv` files in `results/` (by default) and compares them.
  - Default input: `results/` (repo root)
  - Default output: `reports/comparison-predictors/` (repo root)

  Typical way to run:

  ```bash
  python3 scripts/analysis/compare_predictors.py
  ```

  If your system `python3` is missing plotting deps (e.g., `ModuleNotFoundError: matplotlib`), run it using the repo virtualenv:

  ```bash
  ./.venv/bin/python scripts/analysis/compare_predictors.py
  ```

  Useful options:

  ```bash
  # Write outputs somewhere else
  python3 scripts/analysis/compare_predictors.py --output-dir /tmp/cbp_compare

  # Recursively include CSVs in subfolders under results/
  python3 scripts/analysis/compare_predictors.py --recursive

  # Point at a different results folder
  python3 scripts/analysis/compare_predictors.py --results-dir path/to/results
  ```

- `scripts/analysis/generate_delta_figures.py`
  - Expects several predictor CSVs (filenames are hard-coded in `FILES`) to be in the **current working directory**.
  - Writes: `subsumption_analysis/figures/delta_violin_summary_readable.png` (relative to where you run it).

  Example (run from `results/` where the CSVs already exist):

  ```bash
  cd results
  python3 ../scripts/analysis/generate_delta_figures.py
  ```

- `scripts/analysis/delta_curve_visualization.py`
  - Expects a **combined** CSV with columns: `trace_name`, `predictor_name`, `MPKI`.
  - Writes multiple curve plots + per-pair delta CSVs into `--output-dir` (default: `delta_curves/`).

  Example:

  ```bash
  python3 scripts/analysis/delta_curve_visualization.py path/to/combined.csv --output-dir reports/delta_curves
  ```

---

## Quick “where did my PNG go?” checklist

- If you ran `./runPythonScripts.sh baseline`, look under `reports/baseline/**/graphs/`.
- If you ran one of the multi-predictor scripts manually, check the output paths above (many are relative to your current directory).
- If a script says it can’t find `results.csv`, confirm it exists at one of:
  - `./results.csv`
  - `results/results.csv`
  - `scripts/analysis/report_generators/results.csv`

finally: if you’re ever unsure, search under `reports/` for `*.png` and you’ll usually spot the folder you want.
