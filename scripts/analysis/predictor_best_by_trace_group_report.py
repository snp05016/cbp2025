#!/usr/bin/env python3
"""Generate a text report ranking predictors per trace group.

Reads all CSV files in the repo's `results/` directory (or a user-specified directory)
and ranks predictors for each trace group (int/fp/web/infra/compress/media/...).

Ranking metrics:
  - CycWPPKI (lower is better)
  - MPKI     (lower is better)

The report includes per-group aggregated statistics (count/mean/median/std/min/max)
for both metrics, plus a simple "win count" showing how many traces each predictor
is best on (by min CycWPPKI and min MPKI).

Usage:
  python scripts/analysis/predictor_best_by_trace_group_report.py \
      --results-dir results \
      --output reports/predictor_best_by_trace_group.txt

Optional:
  --recursive            Include CSVs under subfolders of results-dir
  --include-per-trace    Emit per-trace winner tables (can be long)

"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


PREDICTOR_LABEL_OVERRIDES: dict[str, str] = {
    "baseline": "Baseline",
    "programming-idiom-predictor-mose": "Programming Idiom Predictor (Mose)",
    "multi-perspective-predictor": "Multi-Perspective Predictor (Jimenez)",
    "tage-sc-l-alberto-ros-results": "TAGE-SC-L (Alberto Ros)",
    "tage-scl-andrez-seznec-results": "TAGE-SCL (Andrez Seznec)",
    "load-value-correlator-man-results": "Load-Value Correlator (Man)",
    "register-value-aware-toru-results": "Register-Value Aware (Toru)",
    "bulleseye-predictor": "Bullseye Predictor",
    "ball-fun": "BALL-FAN",
}

# Canonical groups we expect; we still allow others.
KNOWN_GROUPS = {"int", "fp", "web", "infra", "compress", "media"}


def _repo_root() -> Path:
    # scripts/analysis/<file.py> -> scripts/analysis -> scripts -> repo root
    return Path(__file__).resolve().parents[2]


def _discover_csvs(results_dir: Path, recursive: bool) -> list[Path]:
    if recursive:
        candidates = sorted(p for p in results_dir.rglob("*.csv") if p.is_file())
    else:
        candidates = sorted(p for p in results_dir.glob("*.csv") if p.is_file())

    # Skip common derived artifacts
    skip_names = {
        "combined_predictor_results.csv",
        "pairwise_deltas.csv",
        "subsumption_summary.csv",
    }
    return [p for p in candidates if p.name not in skip_names]


def _label_from_csv_path(path: Path) -> str:
    stem = path.stem
    if stem in PREDICTOR_LABEL_OVERRIDES:
        return PREDICTOR_LABEL_OVERRIDES[stem]

    label = stem.replace("_", " ").replace("-", " ")
    label = re.sub(r"\s+", " ", label).strip()
    return label.title() if label else stem


def _extract_trace_group(row: pd.Series) -> str:
    """Determine trace group robustly.

    Prefer Workload when it's one of the known categories; otherwise fall back
    to parsing from Run (prefix before the first underscore).
    """

    workload = str(row.get("Workload", "")).strip().lower()
    run = str(row.get("Run", "")).strip()

    if workload in KNOWN_GROUPS:
        return workload

    # Many files use Workload='traces' and encode category in Run.
    if run:
        prefix = run.split("_", 1)[0].strip().lower()
        if prefix:
            return prefix

    return "unknown"


def _coerce_numeric(series: pd.Series) -> pd.Series:
    if series.dtype == object:
        # e.g. MR columns may contain percent signs, but MPKI/CycWPPKI should be numeric already
        series = series.astype(str).str.replace(",", "", regex=False)
    return pd.to_numeric(series, errors="coerce")


@dataclass(frozen=True)
class MetricStats:
    count: int
    mean: float
    median: float
    std: float
    min: float
    max: float


def _stats(values: pd.Series) -> MetricStats:
    values = values.dropna()
    if len(values) == 0:
        return MetricStats(0, math.nan, math.nan, math.nan, math.nan, math.nan)

    return MetricStats(
        count=int(values.shape[0]),
        mean=float(values.mean()),
        median=float(values.median()),
        std=float(values.std(ddof=1)) if len(values) > 1 else 0.0,
        min=float(values.min()),
        max=float(values.max()),
    )


def _format_float(v: float, width: int = 10, prec: int = 4) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return f"{'-':>{width}}"
    return f"{v:>{width}.{prec}f}"


def _format_int(v: int, width: int = 6) -> str:
    if v is None:
        return f"{'-':>{width}}"
    return f"{v:>{width}d}"


def _winner_by_mean(table: pd.DataFrame, primary_col: str, secondary_col: str) -> str | None:
    if table.empty:
        return None

    sorted_tbl = table.sort_values([primary_col, secondary_col], ascending=[True, True])
    return str(sorted_tbl.iloc[0]["Predictor"])


def _compute_win_counts(
    df: pd.DataFrame,
    group: str,
    predictors: Iterable[str],
    metric: str,
) -> dict[str, int]:
    """Count, for each predictor, how often it is best (min) per trace in a group."""

    group_df = df[df["TraceGroup"] == group][["Predictor", "Run", metric]].dropna()
    if group_df.empty:
        return {p: 0 for p in predictors}

    wins: dict[str, int] = {p: 0 for p in predictors}

    # For each Run, find predictor(s) with the minimum metric. If tie, split credit.
    for run, run_df in group_df.groupby("Run"):
        min_val = run_df[metric].min()
        best = run_df[run_df[metric] == min_val]["Predictor"].tolist()
        if not best:
            continue
        # Split ties evenly by rounding down (keeps report deterministic and conservative)
        # Alternatively could count all ties as wins, but that inflates totals.
        if len(best) == 1:
            wins[best[0]] += 1

    return wins


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rank predictors per trace group using CycWPPKI and MPKI",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default=str(_repo_root() / "results"),
        help="Directory containing predictor CSV files (default: <repo>/results)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=str(_repo_root() / "reports" / "predictor_best_by_trace_group.txt"),
        help="Output text report path",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively search for CSV files under results-dir",
    )
    parser.add_argument(
        "--include-per-trace",
        action="store_true",
        help="Include per-trace winner tables (can be long)",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_path = Path(args.output)

    csvs = _discover_csvs(results_dir, recursive=args.recursive)
    if not csvs:
        raise SystemExit(f"No CSV files found under: {results_dir}")

    frames: list[pd.DataFrame] = []
    predictor_sources: list[tuple[str, str]] = []

    for csv_path in csvs:
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"Skipping unreadable CSV: {csv_path} ({e})")
            continue

        if "Run" not in df.columns:
            print(f"Skipping CSV missing 'Run': {csv_path.name}")
            continue

        # Filter to passing runs if the column exists
        if "Status" in df.columns:
            df = df[df["Status"].astype(str).str.strip().str.lower() == "pass"].copy()

        predictor = _label_from_csv_path(csv_path)
        predictor_sources.append((predictor, csv_path.name))

        # Ensure required metrics exist
        for col in ("MPKI", "CycWPPKI"):
            if col not in df.columns:
                print(f"Skipping CSV missing '{col}': {csv_path.name}")
                df = None
                break
        if df is None:
            continue

        df = df.copy()
        df["TraceGroup"] = df.apply(_extract_trace_group, axis=1)
        df["Predictor"] = predictor
        df["MPKI"] = _coerce_numeric(df["MPKI"])
        df["CycWPPKI"] = _coerce_numeric(df["CycWPPKI"])

        frames.append(df[["Predictor", "TraceGroup", "Run", "MPKI", "CycWPPKI"]])

    if not frames:
        raise SystemExit("No usable predictor CSVs loaded (missing columns or unreadable files)")

    all_df = pd.concat(frames, ignore_index=True)

    groups = sorted(set(all_df["TraceGroup"].dropna().astype(str)))
    predictors = sorted(set(all_df["Predictor"].dropna().astype(str)))

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write("=" * 100 + "\n")
        f.write("PREDICTOR BEST-BY-TRACE-GROUP REPORT\n")
        f.write("=" * 100 + "\n\n")
        f.write(f"Results directory: {results_dir.resolve()}\n")
        f.write(f"CSV files loaded: {len(csvs)}\n")
        f.write(f"Predictors loaded: {len(predictors)}\n")
        f.write(f"Trace groups found: {', '.join(groups)}\n\n")

        f.write("Predictor sources:\n")
        for pred, src in sorted(predictor_sources, key=lambda t: t[0].lower()):
            f.write(f"  - {pred}: {src}\n")
        f.write("\n")

        # Per-group analysis
        for group in groups:
            group_df = all_df[all_df["TraceGroup"] == group]
            if group_df.empty:
                continue

            # Aggregate stats per predictor
            rows = []
            for pred in predictors:
                pred_df = group_df[group_df["Predictor"] == pred]
                mpki_s = _stats(pred_df["MPKI"])
                cwp_s = _stats(pred_df["CycWPPKI"])
                rows.append(
                    {
                        "Predictor": pred,
                        "N": min(mpki_s.count, cwp_s.count),
                        "MeanCycWPPKI": cwp_s.mean,
                        "MedianCycWPPKI": cwp_s.median,
                        "StdCycWPPKI": cwp_s.std,
                        "MinCycWPPKI": cwp_s.min,
                        "MaxCycWPPKI": cwp_s.max,
                        "MeanMPKI": mpki_s.mean,
                        "MedianMPKI": mpki_s.median,
                        "StdMPKI": mpki_s.std,
                        "MinMPKI": mpki_s.min,
                        "MaxMPKI": mpki_s.max,
                    }
                )

            table = pd.DataFrame(rows)
            table = table[table["N"] > 0].copy()

            # Win counts
            win_cwp = _compute_win_counts(all_df, group, predictors, metric="CycWPPKI")
            win_mpki = _compute_win_counts(all_df, group, predictors, metric="MPKI")
            table["Wins(CycWPPKI)"] = table["Predictor"].map(lambda p: win_cwp.get(p, 0))
            table["Wins(MPKI)"] = table["Predictor"].map(lambda p: win_mpki.get(p, 0))

            # Rankings
            by_cwp = table.sort_values(["MeanCycWPPKI", "MeanMPKI"], ascending=[True, True]).reset_index(drop=True)
            by_mpki = table.sort_values(["MeanMPKI", "MeanCycWPPKI"], ascending=[True, True]).reset_index(drop=True)

            best_cwp = _winner_by_mean(table, primary_col="MeanCycWPPKI", secondary_col="MeanMPKI")
            best_mpki = _winner_by_mean(table, primary_col="MeanMPKI", secondary_col="MeanCycWPPKI")

            f.write("-" * 100 + "\n")
            f.write(f"TRACE GROUP: {group}\n")
            f.write("-" * 100 + "\n")
            f.write(f"Winner by mean CycWPPKI: {best_cwp}\n")
            f.write(f"Winner by mean MPKI:     {best_mpki}\n\n")

            # Full numbers (aggregated) table
            header = (
                f"{'Rank':>4}  {'Predictor':<38}  {'N':>6}  "
                f"{'MeanCWP':>10}  {'MedCWP':>10}  {'StdCWP':>10}  {'MinCWP':>10}  {'MaxCWP':>10}  "
                f"{'MeanMPKI':>10}  {'MedMPKI':>10}  {'StdMPKI':>10}  {'MinMPKI':>10}  {'MaxMPKI':>10}  "
                f"{'WinCWP':>6}  {'WinMPKI':>7}\n"
            )
            f.write("Ranking (primary=Mean CycWPPKI, tiebreak=Mean MPKI):\n")
            f.write(header)

            for idx, row in by_cwp.iterrows():
                f.write(
                    f"{idx+1:>4}  {str(row['Predictor'])[:38]:<38}  "
                    f"{_format_int(int(row['N']))}  "
                    f"{_format_float(row['MeanCycWPPKI'])}  {_format_float(row['MedianCycWPPKI'])}  {_format_float(row['StdCycWPPKI'])}  "
                    f"{_format_float(row['MinCycWPPKI'])}  {_format_float(row['MaxCycWPPKI'])}  "
                    f"{_format_float(row['MeanMPKI'])}  {_format_float(row['MedianMPKI'])}  {_format_float(row['StdMPKI'])}  "
                    f"{_format_float(row['MinMPKI'])}  {_format_float(row['MaxMPKI'])}  "
                    f"{_format_int(int(row['Wins(CycWPPKI)']))}  {_format_int(int(row['Wins(MPKI)']), width=7)}\n"
                )

            f.write("\nRanking (primary=Mean MPKI, tiebreak=Mean CycWPPKI):\n")
            f.write(header)
            for idx, row in by_mpki.iterrows():
                f.write(
                    f"{idx+1:>4}  {str(row['Predictor'])[:38]:<38}  "
                    f"{_format_int(int(row['N']))}  "
                    f"{_format_float(row['MeanCycWPPKI'])}  {_format_float(row['MedianCycWPPKI'])}  {_format_float(row['StdCycWPPKI'])}  "
                    f"{_format_float(row['MinCycWPPKI'])}  {_format_float(row['MaxCycWPPKI'])}  "
                    f"{_format_float(row['MeanMPKI'])}  {_format_float(row['MedianMPKI'])}  {_format_float(row['StdMPKI'])}  "
                    f"{_format_float(row['MinMPKI'])}  {_format_float(row['MaxMPKI'])}  "
                    f"{_format_int(int(row['Wins(CycWPPKI)']))}  {_format_int(int(row['Wins(MPKI)']), width=7)}\n"
                )

            if args.include_per_trace:
                f.write("\nPer-trace winners (min CycWPPKI):\n")
                run_metric = group_df[["Run", "Predictor", "CycWPPKI"]].dropna()
                for run_name, run_df in run_metric.groupby("Run"):
                    best_val = run_df["CycWPPKI"].min()
                    best_preds = run_df[run_df["CycWPPKI"] == best_val]["Predictor"].tolist()
                    best_label = ", ".join(best_preds)
                    f.write(f"  - {run_name}: {best_label} (CycWPPKI={best_val:.4f})\n")

                f.write("\nPer-trace winners (min MPKI):\n")
                run_metric = group_df[["Run", "Predictor", "MPKI"]].dropna()
                for run_name, run_df in run_metric.groupby("Run"):
                    best_val = run_df["MPKI"].min()
                    best_preds = run_df[run_df["MPKI"] == best_val]["Predictor"].tolist()
                    best_label = ", ".join(best_preds)
                    f.write(f"  - {run_name}: {best_label} (MPKI={best_val:.4f})\n")

            f.write("\n")

    print(f"Wrote report: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
