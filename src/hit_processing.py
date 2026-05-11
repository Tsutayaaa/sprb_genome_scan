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
    group_cols = ["target"]
    if "genome_name" in df.columns:
        group_cols = ["genome_name", "target"]
    if "assembly_accession" in df.columns and "assembly_accession" not in group_cols:
        group_cols = ["genome_name", "assembly_accession", "target"] if "genome_name" in df.columns else ["assembly_accession", "target"]
    if "target_fasta" in df.columns and "target_fasta" not in group_cols:
        if "assembly_accession" in df.columns and "genome_name" in df.columns:
            group_cols = ["genome_name", "assembly_accession", "target_fasta", "target"]
        elif "genome_name" in df.columns:
            group_cols = ["genome_name", "target_fasta", "target"]
        else:
            group_cols = ["target_fasta", "target"]

    for group_key, sub in df.groupby(group_cols):
        sub = sub.sort_values("tmin").copy()
        if isinstance(group_key, tuple):
            if len(group_cols) == 4:
                genome_name, assembly_accession, target_fasta, target = group_key
            elif len(group_cols) == 3:
                genome_name, target_fasta, target = group_key
                assembly_accession = sub["assembly_accession"].iloc[0] if "assembly_accession" in sub.columns else ""
            elif len(group_cols) == 2 and group_cols[0] == "genome_name":
                genome_name, target = group_key
                assembly_accession = sub["assembly_accession"].iloc[0] if "assembly_accession" in sub.columns else ""
                target_fasta = sub["target_fasta"].iloc[0] if "target_fasta" in sub.columns else ""
            else:
                target_fasta, target = group_key
                genome_name = sub["genome_name"].iloc[0] if "genome_name" in sub.columns else ""
                assembly_accession = sub["assembly_accession"].iloc[0] if "assembly_accession" in sub.columns else ""
        else:
            target = group_key
            genome_name = sub["genome_name"].iloc[0] if "genome_name" in sub.columns else ""
            assembly_accession = sub["assembly_accession"].iloc[0] if "assembly_accession" in sub.columns else ""
            target_fasta = sub["target_fasta"].iloc[0] if "target_fasta" in sub.columns else ""
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
        row = {
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
        }
        if "genome_name" in df.columns:
            row["genome_name"] = genome_name
        if "assembly_accession" in df.columns:
            row["assembly_accession"] = assembly_accession
        if "organism_name" in df.columns:
            row["organism_name"] = sub["organism_name"].iloc[0]
        if "target_fasta" in df.columns:
            row["target_fasta"] = target_fasta
        rows.append(row)

    out = pd.DataFrame(rows)
    ordered_cols = []
    for col in ["genome_name", "assembly_accession", "organism_name", "target_fasta"]:
        if col in out.columns:
            ordered_cols.append(col)
    ordered_cols.extend([col for col in SUMMARY_COLUMNS if col in out.columns])
    out = out[ordered_cols]
    return out.sort_values(
        [col for col in ["n_unique_families", "union_coverage", "n_hits", "mean_bits"] if col in out.columns],
        ascending=[False, False, False, False][: len([col for col in ["n_unique_families", "union_coverage", "n_hits", "mean_bits"] if col in out.columns])],
    )
