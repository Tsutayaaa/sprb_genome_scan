from __future__ import annotations

from pathlib import Path

import pandas as pd

from .candidate_filter import filter_candidate_proteins
from .config import RunConfig, config_to_dict
from .hit_processing import summarize_by_protein
from .plotting import run_plotting
from .scanner import build_scan_paths, run_scan
from .utils import get_logger


def write_run_info(config: RunConfig, out_file: str | Path) -> None:
    data = config_to_dict(config)
    lines = []
    for key, value in data.items():
        lines.append(f"{key}={value}")
    Path(out_file).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(config: RunConfig) -> dict[str, pd.DataFrame]:
    logger = get_logger(config.verbose)
    filtered_df, plot_df = run_scan(config, logger=logger)
    summary_df = summarize_by_protein(filtered_df)
    summary_with_filter_df = filter_candidate_proteins(summary_df, config.filter_protein)
    candidate_df = summary_with_filter_df[summary_with_filter_df["pass_candidate_filter"]].copy()

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = build_scan_paths(output_dir)
    paths.logs_dir.mkdir(parents=True, exist_ok=True)

    plot_df.to_csv(paths.hits_out, sep="\t", index=False)
    summary_with_filter_df.to_csv(paths.summary_out, sep="\t", index=False)
    candidate_df.to_csv(paths.candidates_out, sep="\t", index=False)
    write_run_info(config, paths.logs_dir / "run_info.txt")

    if config.plot.enabled:
        run_plotting(candidate_df, plot_df, config.plot, output_dir)

    return {
        "hits_df": plot_df,
        "protein_summary": summary_with_filter_df,
        "candidate_df": candidate_df,
    }
