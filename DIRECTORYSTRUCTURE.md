# cbp2025 — directory structure (what’s where)

This repo is a little “pipeline-y”, so the fastest way to orient yourself is: **inputs live in `traces/` and `results/`, scripts live in `scripts/` and the repo root, and outputs land in `reports/`.**

also, if you only read one thing before running stuff, read `instructions.md` first — it explains the intended run flow and where outputs show up.

---

## quick map (top level)

Here’s the main layout, with the folders you’ll actually touch day-to-day.

- `instructions.md`
  - the “how to run everything” doc. it explains the two analysis pipelines (single-predictor vs multi-predictor) and where plots get written.

- `Makefile`
  - build entrypoint for the simulator / predictor interface pieces.

- `cbp` (binary)
  - compiled simulator executable (appears after building).

- `*.sh` helper scripts in repo root
  - `generateReportForPredictor.sh` → generate a report tree for one predictor.
  - `runPythonScripts.sh` → wrapper around single-predictor plotting.
  - `runPredictorAndCompare.sh` → runs predictor report + comparison workflow.
  - `runScripts.sh` → convenience runner (repo-specific).

- `scripts/`
  - python + shell scripts for report generation and comparisons.

- `results/`
  - a collection of per-predictor CSVs (the “inputs” for comparison plots).

- `reports/`
  - the big output tree: graphs, text reports, and compiled HTML artifacts.

- `traces/`
  - trace bundles + documentation, organized by workload category.

- `lib/`
  - the simulation/library code that the predictors plug into.

- `Submissions/`
  - copies of various predictor implementations (and their support headers).

- `analysis/`
  - generated CSVs and derived artifacts for deeper analysis (including subsumption).

- `Assignment2/`
  - LaTeX + compiled report outputs (PDF/HTML) and the figures used in them.

- `papers/`
  - reference papers / writeups that informed some of the included predictors.

---

## the “useful files” (when you’re trying to *do* something)

### if you’re trying to run scripts

Start with `instructions.md`. seriously. it’s the closest thing to a contract for: what to run, where `results.csv` needs to be, and where the PNGs go.

Common entrypoints:

- `generateReportForPredictor.sh`
  - single predictor → generates graphs + text reports under `reports/<predictor_name>/...`.

- `runPredictorAndCompare.sh`
  - one predictor report + then multi-predictor comparisons.

- `scripts/analysis/compare_predictors.py`
  - compares all `*.csv` files under `results/` (by default) and writes to `reports/comparison-predictors/`.

- `scripts/analysis/report_generators/`
  - the individual plotting scripts that get invoked by the single-predictor pipeline.

and yeah, `runPythonScripts.sh` exists too — it’s a wrapper that helps locate/copy a `results.csv` into the right working directory before running the report generators.

### if you’re trying to find traces

There are two trace-ish places:

- `traces/`
  - the “main” distribution folder.
  - contains workload subfolders like `int/`, `fp/`, `web/`, `media/`, `infra/`, plus `.tar.xz` bundles.
  - also includes `readme.txt` and a PDF (`CBP2025-data-dependent-branch-profiles.pdf`).

- `data/traces/`
  - smaller sample traces under `data/traces/sample_traces/` (good for quick sanity runs).

### if you’re trying to find report outputs

Most generated artifacts end up in:

- `reports/`
  - organized by analysis category folders like `01_misprediction_analysis/`, `04_performance_metrics/`, etc.
  - the per-predictor outputs usually look like `reports/<predictor_name>/...`.

A handy “final compiled” snapshot also exists:

- `reports/compiledHtmlReport/`
  - contains `finalCompileReport.html` and a bunch of key PNGs.

There’s also a LaTeX-driven deliverable folder:

- `Assignment2/`
  - includes `finalCompileReport.tex`, `finalCompileReport.pdf`, and the plot images used there.

---

## results & CSVs (inputs to plotting)

You’ll see predictor CSVs in a couple places:

- `results/`
  - contains the comparison-ready CSVs (example: `baseline.csv`, `multi-perspective-predictor.csv`, etc.).
  - also has `baseline/` and `baseline_sim/` subfolders with additional outputs.

- `data/results/`
  - another “results bundle” directory with common predictor CSVs.

- `ReportGenerators/results.csv`
  - a standalone CSV placed with the report-generator tooling (useful as a default/sample).

One detail that trips people up: a lot of the single-predictor plotting scripts expect a file literally named `results.csv` in the current working directory (and `instructions.md` explains how the wrapper script handles that).

---

## source code (.cc / .h) — what’s core vs copied

### core interface + predictor (repo root)

These are the files you’ll usually edit if you’re implementing your own predictor:

- `cond_branch_predictor_interface.cc`
  - interface glue between the simulator and a conditional branch predictor implementation.

- `my_cond_branch_predictor.cc`
- `my_cond_branch_predictor.h`
  - your custom conditional branch predictor implementation.

Core headers at the root (used by predictors / interface code):

- `cbp.h`
- `cbp2016_tage_sc_l.h`
- `cbp2016_tage_sc_l_192kb.h`

### simulator / library code (`lib/`)

the simulator and support code lives here:

- `lib/bp.cc`, `lib/bp.h`
- `lib/cache.cc`, `lib/cache.h`
- `lib/cbp.cc`
- `lib/gzstream.cc`, `lib/gzstream.h`
- `lib/parameters.cc`, `lib/parameters.h`
- `lib/resource_schedule.cc`, `lib/resource_schedule.h`
- `lib/uarchsim.cc`, `lib/uarchsim.h`

and supporting headers (trace readers, interfaces, utilities):

- `lib/trace_reader.h`
- `lib/value_predictor_interface.h`
- `lib/my_value_predictor.cc`, `lib/my_value_predictor.h`
- `lib/stride_prefetcher.h`, `lib/sim_common_structs.h`, `lib/ittage.h`, `lib/fifo.h`, ...

`lib/spdlog/` is a vendored logging library; you normally don’t need to touch it unless you’re debugging prints.

### submissions / third-party predictor drops (`Submissions/`)

This folder contains many predictor implementations (often with their own copy of:
`cbp.h`, `cbp2016_tage_sc_l*.h`, and `my_cond_branch_predictor.h`).

it’s useful for reference and comparison, but it’s also the reason you’ll see a *ton* of duplicated headers.

---

## analysis artifacts

- `analysis/individual/`
  - per-predictor analysis bundles (often zipped), plus extracted working folders.

- `analysis/subsumption/`
  - subsumption/comparison outputs like:
    - `combined_predictor_results.csv`
    - `pairwise_deltas.csv`
    - `subsumption_summary.csv`
    - `figures/`

- `scripts/analysis/`
  - multi-predictor analysis scripts (`compare_predictors.py`, delta tooling, etc.).

---

## small notes (things you’ll notice)

Some folders are intentionally “outputs” and may get large (`reports/`, `analysis/`, and sometimes `results/`).

If you’re lost, a pretty good trick is to just search inside `reports/` for `*.png` — you’ll usually spot the folder the pipeline wrote to.
