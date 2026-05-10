from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_sprb_module_table(
    module_table_path: str | Path,
    cluster_assignments_path: str | Path | None = None,
    cluster_summary_path: str | Path | None = None,
) -> pd.DataFrame:
    module_table_path = Path(module_table_path)
    if not module_table_path.exists():
        raise FileNotFoundError(f"SPRB_MODULE_TABLE not found: {module_table_path}")

    df = pd.read_csv(module_table_path, sep="\t")

    required_std = {"module", "start", "end", "family", "order_index"}
    if required_std.issubset(df.columns):
        out = df.copy()
        out["start"] = pd.to_numeric(out["start"], errors="coerce")
        out["end"] = pd.to_numeric(out["end"], errors="coerce")
        out["order_index"] = pd.to_numeric(out["order_index"], errors="coerce")
        out = out.dropna(subset=["start", "end", "order_index", "family", "module"]).copy()
        return out[["module", "start", "end", "family", "order_index"]]

    required_domain = {"index", "start", "end", "output_label"}
    missing = required_domain - set(df.columns)
    if missing:
        raise ValueError(
            "Input module table is neither standard nor compatible domain table.\n"
            f"Missing columns: {sorted(missing)}"
        )

    if cluster_assignments_path is None or cluster_summary_path is None:
        raise ValueError(
            "cluster_assignments_path and cluster_summary_path are required for domain_table input."
        )

    cluster_assignments_path = Path(cluster_assignments_path)
    cluster_summary_path = Path(cluster_summary_path)
    if not cluster_assignments_path.exists():
        raise FileNotFoundError(f"cluster_assignments.tsv not found: {cluster_assignments_path}")
    if not cluster_summary_path.exists():
        raise FileNotFoundError(f"cluster_summary.tsv not found: {cluster_summary_path}")

    cluster_assign_df = pd.read_csv(cluster_assignments_path, sep="\t")
    cluster_summary_df = pd.read_csv(cluster_summary_path, sep="\t")
    if not {"sequence_id_short", "cluster"}.issubset(cluster_assign_df.columns):
        raise ValueError("cluster_assignments.tsv must contain: sequence_id_short, cluster")
    if not {"cluster_id", "group_name"}.issubset(cluster_summary_df.columns):
        raise ValueError("cluster_summary.tsv must contain: cluster_id, group_name")

    label_to_cluster = dict(zip(cluster_assign_df["sequence_id_short"], cluster_assign_df["cluster"]))
    cluster_to_group = dict(zip(cluster_summary_df["cluster_id"], cluster_summary_df["group_name"]))

    out = df.copy()
    out["module"] = out["output_label"]
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out["order_index"] = pd.to_numeric(out["index"], errors="coerce") - 1
    out["sequence_id_short"] = out["output_label"]
    out["cluster"] = out["sequence_id_short"].map(label_to_cluster)
    out["family"] = out["cluster"].map(cluster_to_group)
    out = out.dropna(subset=["start", "end", "order_index", "module", "family"]).copy()
    return out[["module", "start", "end", "family", "order_index"]]

