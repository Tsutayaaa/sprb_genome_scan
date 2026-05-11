from __future__ import annotations

from collections import Counter

import pandas as pd

from .config import BlockMergeConfig, FilterProteinConfig


SUMMARY_COLUMNS = [
    "target", "tlen", "n_hits", "n_blocks", "n_unique_queries", "n_unique_families", "n_unique_members",
    "cluster_span_start", "cluster_span_end", "cluster_span_len", "cluster_span_coverage",
    "union_covered_len", "union_coverage", "longest_block_len", "max_internal_gap",
    "mean_pident", "mean_bits", "max_bits",
    "best_family_order_string", "best_order_index_string", "linearity_score", "n_order_inversions",
    "repeat_pattern_found", "repeat_motif", "repeat_count", "repeat_motif_length",
    "has_repeated_family", "n_repeated_families", "repeated_family_list",
    "family_list", "member_list", "query_list",
]


BLOCK_COLUMNS = [
    "genome_name", "assembly_accession", "organism_name", "target_fasta",
    "target", "tlen", "block_id", "block_start", "block_end", "block_len",
    "n_hits_in_block", "families_hit", "members_hit", "query_list",
    "best_query", "best_family", "best_member", "best_order_index", "best_order_source",
    "mean_pident", "mean_bits", "max_bits",
]


def _group_cols(df: pd.DataFrame) -> list[str]:
    group_cols = ["target"]
    for col in ["genome_name", "assembly_accession", "organism_name", "target_fasta"]:
        if col in df.columns:
            group_cols.insert(0, col)
    return group_cols


def _metadata_from_group_key(group_key: object, group_cols: list[str], sub: pd.DataFrame) -> dict[str, str]:
    if not isinstance(group_key, tuple):
        group_key = (group_key,)
    key_map = dict(zip(group_cols, group_key))
    return {
        "genome_name": key_map.get("genome_name", sub["genome_name"].iloc[0] if "genome_name" in sub.columns else ""),
        "assembly_accession": key_map.get("assembly_accession", sub["assembly_accession"].iloc[0] if "assembly_accession" in sub.columns else ""),
        "organism_name": key_map.get("organism_name", sub["organism_name"].iloc[0] if "organism_name" in sub.columns else ""),
        "target_fasta": key_map.get("target_fasta", sub["target_fasta"].iloc[0] if "target_fasta" in sub.columns else ""),
        "target": key_map.get("target", sub["target"].iloc[0]),
    }


def _overlap_rate(start_a: int, end_a: int, start_b: int, end_b: int) -> float:
    overlap_start = max(start_a, start_b)
    overlap_end = min(end_a, end_b)
    if overlap_end < overlap_start:
        return 0.0
    overlap_len = overlap_end - overlap_start + 1
    len_a = end_a - start_a + 1
    len_b = end_b - start_b + 1
    return overlap_len / min(len_a, len_b)


def build_reference_lookup(sprb_df: pd.DataFrame) -> dict[str, dict[str, object]]:
    member_to_order: dict[str, int] = {}
    family_to_orders: dict[str, list[int]] = {}

    for _, row in sprb_df.sort_values("order_index").iterrows():
        module = str(row["module"])
        family = str(row["family"])
        order_index = int(row["order_index"])
        family_to_orders.setdefault(family, []).append(order_index)

        aliases = {module}
        if module.startswith("domain_"):
            aliases.add(module.removeprefix("domain_"))
        else:
            aliases.add(f"domain_{module}")
        for alias in aliases:
            member_to_order[alias] = order_index

    family_to_order = {family: min(orders) for family, orders in family_to_orders.items()}
    return {
        "member_to_order": member_to_order,
        "family_to_order": family_to_order,
    }


def _best_hit_reference(hit_row: pd.Series, reference_lookup: dict[str, dict[str, object]] | None) -> tuple[str, str, str, int | None]:
    best_query = str(hit_row["query"])
    best_family = str(hit_row["query_family"])
    best_member = str(hit_row["query_member"])
    if reference_lookup is None:
        return best_query, best_family, best_member, None

    member_to_order = reference_lookup["member_to_order"]
    family_to_order = reference_lookup["family_to_order"]
    if best_member in member_to_order:
        return best_query, best_family, best_member, int(member_to_order[best_member])
    if best_family in family_to_order:
        return best_query, best_family, best_member, int(family_to_order[best_family])
    return best_query, best_family, best_member, None


def merge_hits_to_blocks(
    df: pd.DataFrame,
    block_cfg: BlockMergeConfig,
    reference_lookup: dict[str, dict[str, object]] | None = None,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=BLOCK_COLUMNS)

    rows = []
    grouping_cols = _group_cols(df)

    for group_key, sub in df.groupby(grouping_cols):
        sub = sub.sort_values(["tmin", "tmax", "bits"], ascending=[True, True, False]).copy()
        meta = _metadata_from_group_key(group_key, grouping_cols, sub)
        blocks: list[dict] = []

        for _, hit in sub.iterrows():
            hit_start = int(hit["tmin"])
            hit_end = int(hit["tmax"])
            merged = False
            if blocks:
                last = blocks[-1]
                rate = _overlap_rate(last["block_start"], last["block_end"], hit_start, hit_end)
                if rate >= block_cfg.overlap_rate_threshold:
                    last["block_start"] = min(last["block_start"], hit_start)
                    last["block_end"] = max(last["block_end"], hit_end)
                    last["hits"].append(hit)
                    merged = True
            if not merged:
                blocks.append({"block_start": hit_start, "block_end": hit_end, "hits": [hit]})

        for block_index, block in enumerate(blocks, start=1):
            block_hits = pd.DataFrame(block["hits"]).sort_values(["bits", "pident", "alnlen"], ascending=[False, False, False])
            best_hit = block_hits.iloc[0]
            best_query, best_family, best_member, best_order_index = _best_hit_reference(best_hit, reference_lookup)

            families = block_hits["query_family"].astype(str).tolist()
            members = block_hits["query_member"].astype(str).tolist()
            queries = block_hits["query"].astype(str).tolist()
            row = {
                "target": meta["target"],
                "tlen": int(block_hits["tlen"].max()),
                "block_id": block_index,
                "block_start": block["block_start"],
                "block_end": block["block_end"],
                "block_len": block["block_end"] - block["block_start"] + 1,
                "n_hits_in_block": len(block_hits),
                "families_hit": ",".join(dict.fromkeys(families)),
                "members_hit": ",".join(dict.fromkeys(members)),
                "query_list": ",".join(queries),
                "best_query": best_query,
                "best_family": best_family,
                "best_member": best_member,
                "best_order_index": best_order_index,
                "best_order_source": "member" if best_order_index is not None and best_member in (reference_lookup or {}).get("member_to_order", {}) else ("family" if best_order_index is not None else "missing"),
                "mean_pident": block_hits["pident"].mean(),
                "mean_bits": block_hits["bits"].mean(),
                "max_bits": block_hits["bits"].max(),
            }
            for col in ["genome_name", "assembly_accession", "organism_name", "target_fasta"]:
                if meta[col]:
                    row[col] = meta[col]
            rows.append(row)

    out = pd.DataFrame(rows)
    ordered_cols = [col for col in BLOCK_COLUMNS if col in out.columns]
    return out[ordered_cols].sort_values(
        [col for col in ["genome_name", "assembly_accession", "target", "block_start", "block_end"] if col in out.columns]
    )


def _linearity_metrics(order_values: list[int | None]) -> tuple[float, int]:
    valid_orders = [value for value in order_values if value is not None]
    if len(valid_orders) <= 1:
        return 1.0, 0

    nondecreasing = 0
    inversions = 0
    for prev, curr in zip(valid_orders, valid_orders[1:]):
        if curr >= prev:
            nondecreasing += 1
        else:
            inversions += 1
    score = nondecreasing / (len(valid_orders) - 1)
    return score, inversions


def _detect_repeat_pattern(
    family_sequence: list[str],
    min_len: int,
    max_len: int,
    min_occurrences: int,
    max_gap_blocks: int,
) -> dict[str, object]:
    n = len(family_sequence)
    best_result = {
        "repeat_pattern_found": False,
        "repeat_motif": "",
        "repeat_count": 0,
        "repeat_motif_length": 0,
    }
    if n < max(2, min_len):
        return best_result

    for motif_len in range(min_len, min(max_len, n) + 1):
        for start in range(0, n - motif_len + 1):
            motif = family_sequence[start:start + motif_len]
            if len(motif) < motif_len:
                continue

            positions = [start]
            cursor = start + motif_len
            while cursor <= n - motif_len:
                found = False
                for gap in range(0, max_gap_blocks + 1):
                    candidate_start = cursor + gap
                    if candidate_start > n - motif_len:
                        continue
                    if family_sequence[candidate_start:candidate_start + motif_len] == motif:
                        positions.append(candidate_start)
                        cursor = candidate_start + motif_len
                        found = True
                        break
                if not found:
                    break

            if len(positions) >= min_occurrences:
                result = {
                    "repeat_pattern_found": True,
                    "repeat_motif": ",".join(motif),
                    "repeat_count": len(positions),
                    "repeat_motif_length": motif_len,
                }
                if (
                    result["repeat_count"] > best_result["repeat_count"]
                    or (
                        result["repeat_count"] == best_result["repeat_count"]
                        and result["repeat_motif_length"] > best_result["repeat_motif_length"]
                    )
                ):
                    best_result = result
    return best_result


def summarize_by_protein(blocks_df: pd.DataFrame, filter_cfg: FilterProteinConfig | None = None) -> pd.DataFrame:
    if blocks_df.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    if filter_cfg is None:
        filter_cfg = FilterProteinConfig()

    rows = []
    group_cols = _group_cols(blocks_df)

    for group_key, sub in blocks_df.groupby(group_cols):
        sub = sub.sort_values("block_start").copy()
        meta = _metadata_from_group_key(group_key, group_cols, sub)

        tlen = int(sub["tlen"].max())
        cluster_start = int(sub["block_start"].min())
        cluster_end = int(sub["block_end"].max())
        cluster_len = cluster_end - cluster_start + 1
        cluster_span_coverage = cluster_len / tlen if tlen > 0 else 0.0
        union_covered_len = int(sub["block_len"].sum())
        union_coverage = union_covered_len / tlen if tlen > 0 else 0.0
        n_blocks = len(sub)
        longest_block_len = int(sub["block_len"].max()) if len(sub) > 0 else 0

        gaps = []
        prev_end = None
        for _, row in sub.iterrows():
            if prev_end is not None:
                gaps.append(int(row["block_start"]) - int(prev_end))
            prev_end = int(row["block_end"])
        max_internal_gap = max(gaps) if gaps else 0

        family_tokens = []
        member_tokens = []
        query_tokens = []
        best_family_sequence = []
        best_order_values = []
        for _, row in sub.iterrows():
            family_tokens.extend([x for x in str(row["families_hit"]).split(",") if x])
            member_tokens.extend([x for x in str(row["members_hit"]).split(",") if x])
            query_tokens.extend([x for x in str(row["query_list"]).split(",") if x])
            best_family_sequence.append(str(row["best_family"]))
            order_value = row["best_order_index"]
            if pd.isna(order_value):
                best_order_values.append(None)
            else:
                best_order_values.append(int(order_value))

        fam_counts = Counter(family_tokens)
        repeated_fams = sorted([family for family, count in fam_counts.items() if count >= 2])
        linearity_score, n_order_inversions = _linearity_metrics(best_order_values)
        repeat_info = _detect_repeat_pattern(
            family_sequence=best_family_sequence,
            min_len=filter_cfg.repeat_min_motif_length,
            max_len=filter_cfg.repeat_max_motif_length,
            min_occurrences=filter_cfg.repeat_min_occurrences,
            max_gap_blocks=filter_cfg.repeat_max_gap_blocks,
        )

        row = {
            "target": meta["target"],
            "tlen": tlen,
            "n_hits": int(sub["n_hits_in_block"].sum()),
            "n_blocks": n_blocks,
            "n_unique_queries": len(set(query_tokens)),
            "n_unique_families": len(set(family_tokens)),
            "n_unique_members": len(set(member_tokens)),
            "cluster_span_start": cluster_start,
            "cluster_span_end": cluster_end,
            "cluster_span_len": cluster_len,
            "cluster_span_coverage": cluster_span_coverage,
            "union_covered_len": union_covered_len,
            "union_coverage": union_coverage,
            "longest_block_len": longest_block_len,
            "max_internal_gap": max_internal_gap,
            "mean_pident": sub["mean_pident"].mean(),
            "mean_bits": sub["mean_bits"].mean(),
            "max_bits": sub["max_bits"].max(),
            "best_family_order_string": ",".join(best_family_sequence),
            "best_order_index_string": ",".join("" if value is None else str(value) for value in best_order_values),
            "linearity_score": linearity_score,
            "n_order_inversions": n_order_inversions,
            "repeat_pattern_found": repeat_info["repeat_pattern_found"],
            "repeat_motif": repeat_info["repeat_motif"],
            "repeat_count": repeat_info["repeat_count"],
            "repeat_motif_length": repeat_info["repeat_motif_length"],
            "has_repeated_family": len(repeated_fams) > 0,
            "n_repeated_families": len(repeated_fams),
            "repeated_family_list": ",".join(repeated_fams),
            "family_list": ",".join(sorted(set(family_tokens))),
            "member_list": ",".join(sorted(set(member_tokens))),
            "query_list": ",".join(query_tokens),
        }
        for col in ["genome_name", "assembly_accession", "organism_name", "target_fasta"]:
            if meta[col]:
                row[col] = meta[col]
        rows.append(row)

    out = pd.DataFrame(rows)
    ordered_cols = []
    for col in ["genome_name", "assembly_accession", "organism_name", "target_fasta"]:
        if col in out.columns:
            ordered_cols.append(col)
    ordered_cols.extend([col for col in SUMMARY_COLUMNS if col in out.columns])
    out = out[ordered_cols]
    return out.sort_values(
        ["union_coverage", "linearity_score", "repeat_count", "n_blocks", "mean_bits"],
        ascending=[False, False, False, False, False],
    )
