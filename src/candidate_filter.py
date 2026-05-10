from __future__ import annotations

import pandas as pd

from .config import FilterProteinConfig


def filter_candidate_proteins(summary_df: pd.DataFrame, filter_cfg: FilterProteinConfig) -> pd.DataFrame:
    if summary_df.empty:
        out = summary_df.copy()
        out["pass_candidate_filter"] = pd.Series(dtype=bool)
        return out

    out = summary_df.copy()
    out["pass_candidate_filter"] = (
        (out["tlen"] >= filter_cfg.min_length) &
        (out["n_hits"] >= filter_cfg.min_hits_per_protein) &
        (out["n_unique_families"] >= filter_cfg.min_unique_families) &
        (out["n_unique_members"] >= filter_cfg.min_unique_members) &
        (out["cluster_span_coverage"] >= filter_cfg.min_cluster_span_coverage) &
        (out["union_covered_len"] >= filter_cfg.min_union_covered_len) &
        (out["union_coverage"] >= filter_cfg.min_union_coverage)
    )
    if filter_cfg.require_repeated_family:
        out["pass_candidate_filter"] &= out["has_repeated_family"]
    return out

