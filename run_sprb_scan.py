#!/usr/bin/env python3
from __future__ import annotations

import argparse

from src.config import load_config
from src.pipeline import run_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description="Run modular SprB-like genome scan.")
    parser.add_argument("config", help="Path to JSON or YAML config file.")
    args = parser.parse_args()

    config = load_config(args.config)
    results = run_pipeline(config)

    print(f"[OK] output dir: {config.output_dir}")
    print(f"[OK] results:")
    print(f"  - hits.tsv ({len(results['hits_df'])} rows)")
    print(f"  - protein_hit_summary.tsv ({len(results['protein_summary'])} rows)")
    print(f"  - sprb_like_candidates.tsv ({len(results['candidate_df'])} rows)")
    if config.plot.enabled:
        print("  - plots/")


if __name__ == "__main__":
    main()

