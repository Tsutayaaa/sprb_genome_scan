#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import pandas as pd

from src.config import load_scan_config
from src.candidate_filter import filter_candidate_proteins
from src.hit_processing import summarize_by_protein
from src.scanner import build_scan_paths, discover_targets, run_scan_for_target
from src.utils import get_logger


DEFAULT_SCAN_CONFIG_PATH = Path(__file__).resolve().parent / "configs" / "scan_config.json"


def write_run_info(config: object, out_file: str | Path) -> None:
    from src.config import config_to_dict

    data = config_to_dict(config)
    lines = [f"{key}={value}" for key, value in data.items()]
    Path(out_file).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    config = load_scan_config(DEFAULT_SCAN_CONFIG_PATH)
    logger = get_logger(config.verbose)
    all_filtered = []
    all_plot = []
    targets = discover_targets(config)
    logger.info("Resolved %d target input(s)", len(targets))

    for target_spec in targets:
        logger.info("Scanning target set: %s", target_spec.genome_name)
        filtered_df, plot_df = run_scan_for_target(config, target_spec, logger=logger)
        summary_df = summarize_by_protein(filtered_df)
        summary_with_filter_df = filter_candidate_proteins(summary_df, config.filter_protein)
        candidate_df = summary_with_filter_df[summary_with_filter_df["pass_candidate_filter"]].copy()

        target_output_dir = Path(target_spec.output_dir)
        target_output_dir.mkdir(parents=True, exist_ok=True)
        target_paths = build_scan_paths(target_output_dir)
        target_paths.logs_dir.mkdir(parents=True, exist_ok=True)

        plot_df.to_csv(target_paths.hits_out, sep="\t", index=False)
        summary_with_filter_df.to_csv(target_paths.summary_out, sep="\t", index=False)
        candidate_df.to_csv(target_paths.candidates_out, sep="\t", index=False)
        write_run_info(config, target_paths.logs_dir / "scan_run_info.txt")

        all_filtered.append(filtered_df)
        all_plot.append(plot_df)

    filtered_df = pd.concat(all_filtered, ignore_index=True) if all_filtered else pd.DataFrame()
    plot_df = pd.concat(all_plot, ignore_index=True) if all_plot else pd.DataFrame()
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
    write_run_info(config, paths.logs_dir / "scan_run_info.txt")

    print(f"[OK] scan config: {DEFAULT_SCAN_CONFIG_PATH}")
    print(f"[OK] output dir: {config.output_dir}")
    print(f"[OK] results:")
    print(f"  - hits.tsv ({len(plot_df)} rows)")
    print(f"  - protein_hit_summary.tsv ({len(summary_with_filter_df)} rows)")
    print(f"  - sprb_like_candidates.tsv ({len(candidate_df)} rows)")


if __name__ == "__main__":
    main()
