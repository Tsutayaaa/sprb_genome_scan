from __future__ import annotations

import pandas as pd

from .utils import merge_intervals


SUMMARY_COLUMNS = [
    "target", "tlen", "n_hits", "n_unique_queries", "n_unique_families", "n_unique_members",
    "cluster_span_start", "cluster_span_end", "cluster_span_len", "cluster_span_coverage",
    "union_covered_len", "union_coverage", "n_blocks", "longest_block_len", "max_internal_gap",
    "mean_pident", "mean_bits", "max_bits", "has_repeated_family", "n_repeated_families",
    "repeated_family_list", "family_list", "member_list", "query_list",
]


def summarize_by_protein(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    rows = []
    for target, sub in df.groupby("target"):
        sub = sub.sort_values("tmin").copy()
        tlen = int(sub["tlen"].max())
        cluster_start = int(sub["tmin"].min())
        cluster_end = int(sub["tmax"].max())
        cluster_len = cluster_end - cluster_start + 1
        cluster_span_coverage = cluster_len / tlen if tlen > 0 else 0.0

        merged = merge_intervals((int(row["tmin"]), int(row["tmax"])) for _, row in sub.iterrows())
        union_covered_len = sum(end - start + 1 for start, end in merged)
        union_coverage = union_covered_len / tlen if tlen > 0 else 0.0
        n_blocks = len(merged)
        longest_block_len = max((end - start + 1 for start, end in merged), default=0)

        gaps = []
        prev_end = None
        for _, row in sub.iterrows():
            if prev_end is not None:
                gaps.append(int(row["tmin"]) - int(prev_end))
            prev_end = int(row["tmax"])
        max_internal_gap = max(gaps) if gaps else 0

        fam_counts = sub["query_family"].value_counts()
        repeated_fams = sorted(fam_counts[fam_counts >= 2].index.tolist())
        rows.append({
            "target": target,
            "tlen": tlen,
            "n_hits": len(sub),
            "n_unique_queries": sub["query"].nunique(),
            "n_unique_families": sub["query_family"].nunique(),
            "n_unique_members": sub["query_member"].nunique(),
            "cluster_span_start": cluster_start,
            "cluster_span_end": cluster_end,
            "cluster_span_len": cluster_len,
            "cluster_span_coverage": cluster_span_coverage,
            "union_covered_len": union_covered_len,
            "union_coverage": union_coverage,
            "n_blocks": n_blocks,
            "longest_block_len": longest_block_len,
            "max_internal_gap": max_internal_gap,
            "mean_pident": sub["pident"].mean(),
            "mean_bits": sub["bits"].mean(),
            "max_bits": sub["bits"].max(),
            "has_repeated_family": len(repeated_fams) > 0,
            "n_repeated_families": len(repeated_fams),
            "repeated_family_list": ",".join(repeated_fams),
            "family_list": ",".join(sorted(sub["query_family"].unique())),
            "member_list": ",".join(sorted(sub["query_member"].unique())),
            "query_list": ",".join(sub["query"].tolist()),
        })

    out = pd.DataFrame(rows)
    return out.sort_values(
        ["n_unique_families", "union_coverage", "n_hits", "mean_bits"],
        ascending=[False, False, False, False],
    )

