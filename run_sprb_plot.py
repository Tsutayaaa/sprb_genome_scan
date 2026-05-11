#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import pandas as pd

from src.config import load_plot_config
from src.plotting import run_plotting_for_result_dir
from src.utils import get_logger, progress_stats
import time


DEFAULT_PLOT_CONFIG_PATH = Path(__file__).resolve().parent / "configs" / "plot_config.json"


def write_run_info(config: object, out_file: str | Path) -> None:
    from src.config import config_to_dict

    data = config_to_dict(config)
    lines = [f"{key}={value}" for key, value in data.items()]
    Path(out_file).write_text("\n".join(lines) + "\n", encoding="utf-8")


def discover_result_dirs(results_root: Path, results_mode: str) -> list[Path]:
    mode = results_mode.lower()
    if mode not in {"auto", "single_result_dir", "per_genome_dirs"}:
        raise ValueError(f"Unsupported results_mode: {results_mode}")

    required = {"blocks.tsv", "protein_hit_summary.tsv", "sprb_like_candidates.tsv"}
    if all((results_root / name).exists() for name in required):
        if mode in {"auto", "single_result_dir"}:
            return [results_root]

    result_dirs = []
    for child in sorted(results_root.iterdir()):
        if not child.is_dir():
            continue
        if all((child / name).exists() for name in required):
            result_dirs.append(child)

    if result_dirs:
        return result_dirs
    raise FileNotFoundError(f"No result directories containing {sorted(required)} were found under: {results_root}")


def main() -> None:
    config = load_plot_config(DEFAULT_PLOT_CONFIG_PATH)
    logger = get_logger(config.verbose)
    results_root = Path(config.results_root)
    if not results_root.exists():
        raise FileNotFoundError(f"Results root not found: {results_root}")

    result_dirs = discover_result_dirs(results_root, config.results_mode)
    logger.info("Resolved %d result directory(ies)", len(result_dirs))
    run_start = time.time()

    plotted_targets = []
    for index, result_dir in enumerate(result_dirs, start=1):
        before_stats = progress_stats(run_start, index - 1, len(result_dirs))
        logger.info(
            "Progress %d/%d (%.1f%%) | elapsed=%s | eta=%s | plotting=%s",
            index,
            len(result_dirs),
            before_stats["progress_pct"],
            before_stats["elapsed_text"],
            before_stats["eta_text"],
            result_dir.name,
        )
        plotted = run_plotting_for_result_dir(result_dir, Path(config.output_dir), config.plot)
        plotted_targets.extend(plotted)
        after_stats = progress_stats(run_start, index, len(result_dirs))
        logger.info(
            "Completed %d/%d | result_dir=%s | plotted_targets=%d | elapsed=%s | avg/dir=%s | eta=%s",
            index,
            len(result_dirs),
            result_dir.name,
            len(plotted),
            after_stats["elapsed_text"],
            after_stats["avg_text"],
            after_stats["eta_text"],
        )

    logs_dir = Path(config.output_dir) / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    write_run_info(config, logs_dir / "plot_run_info.txt")
    plotted_targets_df = pd.DataFrame({"target": plotted_targets})
    total_stats = progress_stats(run_start, len(result_dirs), len(result_dirs))

    print(f"[OK] plot config: {DEFAULT_PLOT_CONFIG_PATH}")
    print(f"[OK] results root: {config.results_root}")
    print(f"[OK] plot output dir: {config.output_dir}")
    print(f"[OK] elapsed: {total_stats['elapsed_text']}")
    print(f"[OK] plotted targets: {len(plotted_targets_df)}")


if __name__ == "__main__":
    main()
