from __future__ import annotations

import pandas as pd

from .config import FilterProteinConfig


def filter_candidate_proteins(summary_df: pd.DataFrame, filter_cfg: FilterProteinConfig) -> pd.DataFrame:
    if summary_df.empty:
        out = summary_df.copy()
        out["pass_candidate_filter"] = pd.Series(dtype=bool)
        out["failed_reasons"] = pd.Series(dtype=str)
        return out

    out = summary_df.copy()
    out["pass_min_length"] = out["tlen"] >= filter_cfg.min_length
    out["pass_min_blocks"] = out["n_blocks"] >= filter_cfg.min_blocks
    out["pass_unique_families"] = out["n_unique_families"] >= filter_cfg.min_unique_families
    out["pass_unique_members"] = out["n_unique_members"] >= filter_cfg.min_unique_members
    out["pass_cluster_span_coverage"] = out["cluster_span_coverage"] >= filter_cfg.min_cluster_span_coverage
    out["pass_union_covered_len"] = out["union_covered_len"] >= filter_cfg.min_union_covered_len
    out["pass_union_coverage"] = out["union_coverage"] >= filter_cfg.min_union_coverage
    out["pass_linearity_score"] = out["linearity_score"] >= filter_cfg.min_linearity_score
    out["pass_order_inversions"] = out["n_order_inversions"] <= filter_cfg.max_order_inversions
    out["pass_candidate_filter"] = (
        out["pass_min_length"] &
        out["pass_min_blocks"] &
        out["pass_unique_families"] &
        out["pass_unique_members"] &
        out["pass_cluster_span_coverage"] &
        out["pass_union_covered_len"] &
        out["pass_union_coverage"] &
        out["pass_linearity_score"] &
        out["pass_order_inversions"]
    )
    if filter_cfg.require_repeat_pattern:
        out["pass_repeat_pattern"] = out["repeat_pattern_found"]
        out["pass_repeat_count"] = out["repeat_count"] >= filter_cfg.min_repeat_count
        out["pass_candidate_filter"] &= out["pass_repeat_pattern"]
        out["pass_candidate_filter"] &= out["pass_repeat_count"]
    else:
        out["pass_repeat_pattern"] = True
        out["pass_repeat_count"] = True

    check_columns = [
        ("pass_min_length", "min_length"),
        ("pass_min_blocks", "min_blocks"),
        ("pass_unique_families", "unique_families"),
        ("pass_unique_members", "unique_members"),
        ("pass_cluster_span_coverage", "cluster_span_coverage"),
        ("pass_union_covered_len", "union_covered_len"),
        ("pass_union_coverage", "union_coverage"),
        ("pass_linearity_score", "linearity_score"),
        ("pass_order_inversions", "order_inversions"),
        ("pass_repeat_pattern", "repeat_pattern"),
        ("pass_repeat_count", "repeat_count"),
    ]

    def _failed_reasons(row: pd.Series) -> str:
        reasons = [label for column, label in check_columns if column in row.index and not bool(row[column])]
        return ",".join(reasons)

    out["failed_reasons"] = out.apply(_failed_reasons, axis=1)
    return out
