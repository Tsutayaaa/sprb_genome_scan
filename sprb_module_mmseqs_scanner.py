#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
sprb_module_mmseqs_scanner.py

功能
----
使用 SprB 模块 panel fasta，对指定目标蛋白组执行 MMseqs2 检索，
并输出命中结果与按蛋白汇总的候选注释。

绘图后端
--------
使用 pyGenomeViz 绘制类似 clinker 风格的 SprB vs candidate architecture 图：

- 上方：SprB 模块参考轨道
- 下方：目标蛋白命中轨道
- 中间：family 连接线
- 支持 one-to-many / best-link 两种连线模式

兼容输入
--------
1. 标准表：
   module / start / end / family / order_index

2. domain_table.tsv：
   index / start / end / output_label / ...
   再结合：
   - cluster_assignments.tsv
   - cluster_summary.tsv
   自动补 family
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List
import re

import pandas as pd
from pygenomeviz import GenomeViz

# =========================================================
# 直接在这里改参数
# =========================================================

# query panel fasta
QUERY_FASTA = "/Users/shulei/Documents/2026_Chemotaxis/sprb_domains/grouped_for_hmm_blosum62/sprb_smf_panel.fasta"

# 目标蛋白组 fasta
TARGET_FASTA = "/Users/shulei/PycharmProjects/Biopython/sprB/ncbi_gliding_genomes/Flavobacterium_johnsoniae_ATCC_17061/GCF_000016645.1_ASM1664v1_protein.faa"

# 输出目录
OUTPUT_DIR = "/Users/shulei/Documents/2026_Chemotaxis/sprb_domains/mmseqs_scan_example_FJ_2"

# mmseqs 可执行文件
MMSEQS_BIN = "mmseqs"

# 线程数
THREADS = 8

# 是否强制重建数据库和重新搜索
FORCE_RERUN = True

# mmseqs 搜索参数
# -s: 搜索灵敏度，重复模块扫描推荐 5.7-7.5 起步
SEARCH_SENSITIVITY = 7.5
# --min-seq-id: 允许较远的重复单元同源
MIN_SEQ_ID = 0.15
# -e: 显著性阈值，先保持平衡设置
EVALUE = 1e-3
# --max-seqs: 单 query 保留的最大结果数
MAX_SEQS = 5000
# coverage 先保持宽松，避免局部重复位点被提前过滤
COV_MODE = 0
MIN_COVERAGE = 0.0
# 允许同一 query-target 返回多个替代局部对齐
ALT_ALI = 20

# 过滤阈值
MIN_BITS = 30.0
MIN_ALN_LEN = 50
MIN_PIDENT = 15.0

# 候选蛋白判定参数（初步）
MIN_HITS_PER_PROTEIN = 3
MIN_UNIQUE_FAMILIES = 2
MAX_GAP_BETWEEN_HITS = 120

# =========================================================
# 蛋白级过滤参数
# =========================================================

# 候选 target 蛋白的最小长度；过短蛋白通常只是残余同源片段
MIN_TARGET_LENGTH = 500

# 命中总跨度 [min(tmin), max(tmax)] 占 target 全长的最小比例
MIN_CLUSTER_SPAN_COVERAGE = 0.20

# 合并所有 hit 区间后的 union 覆盖率下限，避免只有少量零散局部命中
MIN_UNION_COVERAGE = 0.08

# 合并所有 hit 区间后的净被命中区域总长度下限
MIN_UNION_COVERED_LEN = 200

# target 上至少命中的不同 query member 数
MIN_UNIQUE_MEMBERS = 2

# 是否要求至少存在一个重复出现的 family 才通过蛋白级过滤
REQUIRE_REPEATED_FAMILY = False

# 是否打印详细日志
VERBOSE = True

# =========================================================
# 绘图参数
# =========================================================

ENABLE_PLOTTING = True

# SprB 模块参考表：可直接用 domain_table.tsv
SPRB_MODULE_TABLE = "/Users/shulei/Documents/2026_Chemotaxis/sprb_domains/domain_table.tsv"

# 若不是标准表，则需下面两个文件做 family 映射
CLUSTER_ASSIGNMENTS_TSV = "/Users/shulei/Documents/2026_Chemotaxis/sprb_domains/grouped_for_hmm_blosum62/cluster_assignments.tsv"
CLUSTER_SUMMARY_TSV = "/Users/shulei/Documents/2026_Chemotaxis/sprb_domains/grouped_for_hmm_blosum62/cluster_summary.tsv"

# SprB 全长长度
SPRB_LENGTH = 6497

# 画哪些蛋白
# "top_candidates" / "specified"
PLOT_TARGET_MODE = "top_candidates"

# 如果画 top candidates，画前多少个
TOP_N_CANDIDATES = 10

# 如果用 specified，填 protein id
SPECIFIED_TARGETS = [
    # "WP_011962197.1",
]

# 连线模式
# "all"  : 同 family 的 SprB 模块全部连到目标 hit
# "best" : 优先连到 query 对应的原始模块；找不到则退回 family 的第一个模块
PLOT_LINK_MODE = "all"

# 图像输出格式
PLOT_SAVE_PNG = True
PLOT_SAVE_PDF = True
PLOT_SAVE_HTML = False

# pyGenomeViz 画布参数（v1.6.1）
PLOT_FIG_WIDTH = 18
PLOT_FIG_TRACK_HEIGHT = 1.2
# "left", "center", "right"
PLOT_TRACK_ALIGN_TYPE = "center"
PLOT_FEATURE_TRACK_RATIO = 0.32
PLOT_LINK_TRACK_RATIO = 0.90
PLOT_THEME = "light"
PLOT_SHOW_AXIS = False

# track 样式
TRACK_LABEL_SIZE = 16
TRACK_LABEL_MARGIN = 0.02
TRACK_ALIGN_LABEL = True
TRACK_LINE_KWS = dict(color="black", lw=0.8)

# feature 样式
DRAW_LABELS = True

# "bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"
FEATURE_PLOTSTYLE = "box"
FEATURE_LINEWIDTH = 0.6
FEATURE_LABELSIZE = 8
FEATURE_TEXT_ROTATION = 90

# target 轨道是否显示 family 标签
SHOW_TARGET_LABELS = False

# link 样式
LINK_CURVE = True
LINK_SIZE = 0.9

# 连线透明度映射到相似性
LINK_ALPHA_MIN = 0.12
LINK_ALPHA_MAX = 0.80

# 输出 dpi
PLOT_DPI = 300

# =========================================================
# 绘图筛选参数（只影响显示，不影响主结果）
# =========================================================

PLOT_MIN_BITS = 35.0
PLOT_MIN_ALN_LEN = 60
PLOT_MIN_PIDENT = 15.0

# =========================================================
# 工具函数
# =========================================================

def log(msg: str) -> None:
    if VERBOSE:
        print(msg)

def run_command(cmd: List[str]) -> None:
    log("[CMD] " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(cmd)}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )

def remove_mmseqs_artifacts(prefix: Path) -> None:
    """删除某个 mmseqs 前缀对应的所有文件"""
    parent = prefix.parent
    stem = prefix.name
    for p in parent.glob(stem + "*"):
        if p.is_file():
            p.unlink()

def prepare_for_rerun(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    tmp_dir: Path,
    raw_out: Path,
    filtered_out: Path,
    plot_filtered_out: Path,
    summary_out: Path,
    summary_filtered_out: Path,
    candidate_out: Path,
    family_summary_out: Path,
    by_target_dir: Path,
) -> None:
    """FORCE_RERUN=True 时，先清理旧文件"""
    if not FORCE_RERUN:
        return

    log("[INFO] FORCE_RERUN=True, removing old MMseqs artifacts")

    remove_mmseqs_artifacts(query_db)
    remove_mmseqs_artifacts(target_db)
    remove_mmseqs_artifacts(result_db)

    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    for f in [raw_out, filtered_out, plot_filtered_out, summary_out, summary_filtered_out, candidate_out, family_summary_out]:
        if f.exists():
            f.unlink()
    if by_target_dir.exists():
        shutil.rmtree(by_target_dir)

def ensure_mmseqs_db(fasta: Path, db_path: Path) -> None:
    db_type_file = db_path.with_suffix(".dbtype")
    if db_type_file.exists() and not FORCE_RERUN:
        log(f"[INFO] MMseqs db already exists: {db_path}")
        return

    cmd = [MMSEQS_BIN, "createdb", str(fasta), str(db_path)]
    run_command(cmd)

def run_mmseqs_search(query_db: Path, target_db: Path, result_db: Path, tmp_dir: Path) -> None:
    db_type_file = result_db.with_suffix(".dbtype")
    if db_type_file.exists() and not FORCE_RERUN:
        log(f"[INFO] MMseqs result db already exists: {result_db}")
        return

    cmd = [
        MMSEQS_BIN, "search",
        str(query_db),
        str(target_db),
        str(result_db),
        str(tmp_dir),
        "--threads", str(THREADS),
        "-s", str(SEARCH_SENSITIVITY),
        "-e", str(EVALUE),
        "--min-seq-id", str(MIN_SEQ_ID),
        "--max-seqs", str(MAX_SEQS),
        "--cov-mode", str(COV_MODE),
        "-c", str(MIN_COVERAGE),
        "--alt-ali", str(ALT_ALI),
    ]
    run_command(cmd)

def export_mmseqs_tsv(query_db: Path, target_db: Path, result_db: Path, out_tsv: Path) -> None:
    if out_tsv.exists() and (not FORCE_RERUN):
        log(f"[INFO] MMseqs TSV already exists: {out_tsv}")
        return

    outfmt = ",".join([
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits",
        "qlen", "tlen"
    ])

    cmd = [
        MMSEQS_BIN, "convertalis",
        str(query_db),
        str(target_db),
        str(result_db),
        str(out_tsv),
        "--format-output", outfmt
    ]
    run_command(cmd)

def extract_family(query_id: str) -> str:
    if "__" in query_id:
        return query_id.split("__", 1)[0]
    return query_id

def extract_member(query_id: str) -> str:
    if "__" in query_id:
        return query_id.split("__", 1)[1]
    return query_id

def load_mmseqs_table(path: Path) -> pd.DataFrame:
    cols = [
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits",
        "qlen", "tlen"
    ]
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=cols)

    df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str)

    if len(df) > 0 and str(df.iloc[0]["query"]) == "query":
        df = df.iloc[1:].copy()

    numeric_cols = [
        "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend",
        "evalue", "bits", "qlen", "tlen"
    ]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["fident", "alnlen", "bits", "tstart", "tend", "tlen"]).copy()

    df["pident"] = df["fident"] * 100.0
    df["tmin"] = df[["tstart", "tend"]].min(axis=1)
    df["tmax"] = df[["tstart", "tend"]].max(axis=1)
    df["strand"] = df.apply(lambda r: "+" if r["tend"] >= r["tstart"] else "-", axis=1)
    df["query_family"] = df["query"].map(extract_family)
    df["query_member"] = df["query"].map(extract_member)

    return df

def filter_hits(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    out = df.copy()
    out = out[
        (out["pident"] >= MIN_PIDENT) &
        (out["alnlen"] >= MIN_ALN_LEN) &
        (out["bits"] >= MIN_BITS)
    ].copy()

    return out

def build_plot_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    仅用于绘图的显示筛选，不影响主输出 filtered_df
    """
    if df.empty:
        return df.copy()

    return df[
        (df["bits"] >= PLOT_MIN_BITS) &
        (df["alnlen"] >= PLOT_MIN_ALN_LEN) &
        (df["pident"] >= PLOT_MIN_PIDENT)
    ].copy()

def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []

    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged

def summarize_by_protein(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "target", "tlen", "n_hits", "n_unique_queries", "n_unique_families", "n_unique_members",
            "cluster_span_start", "cluster_span_end", "cluster_span_len", "cluster_span_coverage",
            "union_covered_len", "union_coverage",
            "n_blocks", "longest_block_len", "max_internal_gap",
            "mean_pident", "mean_bits", "max_bits",
            "has_repeated_family", "n_repeated_families", "repeated_family_list",
            "family_list", "member_list", "query_list"
        ])

    rows = []
    for target, sub in df.groupby("target"):
        sub = sub.sort_values("tmin").copy()
        tlen = int(sub["tlen"].max())

        cluster_start = int(sub["tmin"].min())
        cluster_end = int(sub["tmax"].max())
        cluster_len = cluster_end - cluster_start + 1
        cluster_span_coverage = cluster_len / tlen if tlen > 0 else 0.0

        intervals = [(int(r["tmin"]), int(r["tmax"])) for _, r in sub.iterrows()]
        merged = merge_intervals(intervals)
        union_covered_len = sum(e - s + 1 for s, e in merged)
        union_coverage = union_covered_len / tlen if tlen > 0 else 0.0
        n_blocks = len(merged)
        longest_block_len = max((e - s + 1 for s, e in merged), default=0)

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
    out = out.sort_values(
        ["n_unique_families", "union_coverage", "n_hits", "mean_bits"],
        ascending=[False, False, False, False]
    )
    return out

def filter_candidate_proteins(summary_df: pd.DataFrame) -> pd.DataFrame:
    if summary_df.empty:
        out = summary_df.copy()
        out["pass_candidate_filter"] = pd.Series(dtype=bool)
        return out

    out = summary_df.copy()
    out["pass_candidate_filter"] = (
        (out["tlen"] >= MIN_TARGET_LENGTH) &
        (out["n_hits"] >= MIN_HITS_PER_PROTEIN) &
        (out["n_unique_families"] >= MIN_UNIQUE_FAMILIES) &
        (out["n_unique_members"] >= MIN_UNIQUE_MEMBERS) &
        (out["cluster_span_coverage"] >= MIN_CLUSTER_SPAN_COVERAGE) &
        (out["union_covered_len"] >= MIN_UNION_COVERED_LEN) &
        (out["union_coverage"] >= MIN_UNION_COVERAGE)
    )

    if REQUIRE_REPEATED_FAMILY:
        out["pass_candidate_filter"] &= out["has_repeated_family"]

    return out

def summarize_by_family(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "query_family", "n_hits", "n_targets", "mean_pident", "mean_bits"
        ])

    rows = []
    for fam, sub in df.groupby("query_family"):
        rows.append({
            "query_family": fam,
            "n_hits": len(sub),
            "n_targets": sub["target"].nunique(),
            "mean_pident": sub["pident"].mean(),
            "mean_bits": sub["bits"].mean(),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(["n_hits", "n_targets"], ascending=[False, False])
    return out

def detect_clustered_hits(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "target", "n_hits", "n_unique_queries", "n_unique_families",
            "clustered", "max_internal_gap",
            "cluster_span_start", "cluster_span_end", "cluster_span_len",
            "mean_pident", "mean_bits", "family_list", "query_list"
        ])

    rows = []
    for target, sub in df.groupby("target"):
        sub = sub.sort_values("tmin").copy()

        gaps = []
        prev_end = None
        for _, row in sub.iterrows():
            if prev_end is not None:
                gaps.append(int(row["tmin"]) - int(prev_end))
            prev_end = row["tmax"]

        max_gap = max(gaps) if gaps else 0

        n_hits = len(sub)
        n_queries = sub["query"].nunique()
        n_families = sub["query_family"].nunique()
        cluster_start = int(sub["tmin"].min())
        cluster_end = int(sub["tmax"].max())
        cluster_len = cluster_end - cluster_start + 1

        clustered = (
            (n_hits >= MIN_HITS_PER_PROTEIN) and
            (n_families >= MIN_UNIQUE_FAMILIES) and
            (max_gap <= MAX_GAP_BETWEEN_HITS)
        )

        rows.append({
            "target": target,
            "n_hits": n_hits,
            "n_unique_queries": n_queries,
            "n_unique_families": n_families,
            "clustered": clustered,
            "max_internal_gap": max_gap,
            "cluster_span_start": cluster_start,
            "cluster_span_end": cluster_end,
            "cluster_span_len": cluster_len,
            "mean_pident": sub["pident"].mean(),
            "mean_bits": sub["bits"].mean(),
            "family_list": ",".join(sorted(sub["query_family"].unique())),
            "query_list": ",".join(sub["query"].tolist()),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["clustered", "n_hits", "n_unique_families", "mean_bits"],
        ascending=[False, False, False, False]
    )
    return out

def write_run_info(out_file: Path) -> None:
    text = f"""QUERY_FASTA={QUERY_FASTA}
TARGET_FASTA={TARGET_FASTA}
SEARCH_SENSITIVITY={SEARCH_SENSITIVITY}
EVALUE={EVALUE}
MIN_SEQ_ID={MIN_SEQ_ID}
MAX_SEQS={MAX_SEQS}
COV_MODE={COV_MODE}
MIN_COVERAGE={MIN_COVERAGE}
ALT_ALI={ALT_ALI}
MIN_BITS={MIN_BITS}
MIN_ALN_LEN={MIN_ALN_LEN}
MIN_PIDENT={MIN_PIDENT}
PLOT_MIN_BITS={PLOT_MIN_BITS}
PLOT_MIN_ALN_LEN={PLOT_MIN_ALN_LEN}
PLOT_MIN_PIDENT={PLOT_MIN_PIDENT}
MIN_HITS_PER_PROTEIN={MIN_HITS_PER_PROTEIN}
MIN_UNIQUE_FAMILIES={MIN_UNIQUE_FAMILIES}
MAX_GAP_BETWEEN_HITS={MAX_GAP_BETWEEN_HITS}
MIN_TARGET_LENGTH={MIN_TARGET_LENGTH}
MIN_CLUSTER_SPAN_COVERAGE={MIN_CLUSTER_SPAN_COVERAGE}
MIN_UNION_COVERED_LEN={MIN_UNION_COVERED_LEN}
MIN_UNION_COVERAGE={MIN_UNION_COVERAGE}
MIN_UNIQUE_MEMBERS={MIN_UNIQUE_MEMBERS}
REQUIRE_REPEATED_FAMILY={REQUIRE_REPEATED_FAMILY}
THREADS={THREADS}
"""
    out_file.write_text(text, encoding="utf-8")

def write_target_reports(df: pd.DataFrame, summary_df: pd.DataFrame, output_dir: Path) -> None:
    by_target_dir = output_dir / "by_target"
    by_target_dir.mkdir(parents=True, exist_ok=True)

    summary_map = {row["target"]: row for _, row in summary_df.iterrows()}

    for target, sub in df.groupby("target"):
        target_dir = by_target_dir / str(target)
        target_dir.mkdir(parents=True, exist_ok=True)

        sub = sub.sort_values(["tmin", "tmax"]).copy()
        sub.to_csv(target_dir / "hits.tsv", sep="\t", index=False)

        summary_row = summary_map.get(target)
        if summary_row is not None:
            pd.DataFrame([summary_row]).to_csv(target_dir / "summary.tsv", sep="\t", index=False)

        intervals = [(int(r["tmin"]), int(r["tmax"])) for _, r in sub.iterrows()]
        merged = merge_intervals(intervals)
        merged_df = pd.DataFrame([
            {"block_id": i + 1, "start": s, "end": e, "length": e - s + 1}
            for i, (s, e) in enumerate(merged)
        ])
        merged_df.to_csv(target_dir / "merged_blocks.tsv", sep="\t", index=False)

        fam_rows = []
        for fam, fam_sub in sub.groupby("query_family"):
            fam_rows.append({
                "query_family": fam,
                "n_hits": len(fam_sub),
                "mean_bits": fam_sub["bits"].mean(),
                "mean_pident": fam_sub["pident"].mean(),
                "start_min": int(fam_sub["tmin"].min()),
                "end_max": int(fam_sub["tmax"].max()),
                "members": ",".join(sorted(fam_sub["query_member"].unique())),
            })
        fam_df = pd.DataFrame(fam_rows).sort_values(["n_hits", "mean_bits"], ascending=[False, False])
        fam_df.to_csv(target_dir / "family_counts.tsv", sep="\t", index=False)

# =========================================================
# 参考表兼容读取
# =========================================================

def load_sprb_module_table_compatible(
    module_table_path: Path,
    cluster_assignments_path: Path,
    cluster_summary_path: Path
) -> pd.DataFrame:
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
        return out

    required_domain = {"index", "start", "end", "output_label"}
    missing = required_domain - set(df.columns)
    if missing:
        raise ValueError(
            f"Input module table is neither standard nor compatible domain table.\n"
            f"Missing columns: {sorted(missing)}"
        )

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
    out = out[["module", "start", "end", "family", "order_index"]].copy()

    return out

# =========================================================
# pyGenomeViz 绘图函数
# =========================================================

def build_family_color_map(families: list[str]) -> dict[str, str]:
    unique_fams = sorted(set(families))
    base_colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1",
        "#9c755f", "#ff9da7", "#bab0ab", "#86bc86", "#499894",
    ]
    return {fam: base_colors[i % len(base_colors)] for i, fam in enumerate(unique_fams)}

def short_module_name(module_name: str) -> str:
    return re.sub(r"^domain_", "", str(module_name))

def build_ref_rows_for_hit(sprb_df: pd.DataFrame, hit_row: pd.Series) -> pd.DataFrame:
    fam = hit_row["query_family"]
    qmember = hit_row["query_member"]

    if PLOT_LINK_MODE == "all":
        ref_rows = sprb_df[sprb_df["family"] == fam]

    elif PLOT_LINK_MODE == "best":
        ref_rows = sprb_df[
            (sprb_df["family"].astype(str) == str(fam)) &
            (sprb_df["module"].astype(str) == str(qmember))
        ]
        if ref_rows.empty:
            qmember_domain = f"domain_{qmember}" if not str(qmember).startswith("domain_") else str(qmember)
            ref_rows = sprb_df[
                (sprb_df["family"].astype(str) == str(fam)) &
                (sprb_df["module"].astype(str) == qmember_domain)
            ]
        if ref_rows.empty:
            ref_rows = sprb_df[sprb_df["family"] == fam].head(1)
    else:
        raise ValueError(f"Unsupported PLOT_LINK_MODE: {PLOT_LINK_MODE}")

    return ref_rows

def map_similarity_to_alpha(
    value: float,
    vmin: float,
    vmax: float,
    alpha_min: float = LINK_ALPHA_MIN,
    alpha_max: float = LINK_ALPHA_MAX,
) -> float:
    """
    将相似性数值线性映射到透明度
    相似性越高，线越深
    """
    if pd.isna(value):
        return alpha_min
    if vmax <= vmin:
        return alpha_max

    x = (value - vmin) / (vmax - vmin)
    x = max(0.0, min(1.0, x))
    return alpha_min + x * (alpha_max - alpha_min)

def plot_candidate_architecture_pygenomeviz(
    sprb_df: pd.DataFrame,
    hits_df: pd.DataFrame,
    target_id: str,
    outdir: Path,
    output_prefix: str = "sprb_candidate_architecture"
) -> None:
    sub = hits_df[hits_df["target"] == target_id].copy()
    if sub.empty:
        return

    sub = sub.sort_values("tmin").copy()
    target_len = int(sub["tlen"].max())

    all_families = list(sprb_df["family"]) + list(sub["query_family"])
    color_map = build_family_color_map(all_families)

    # 当前 target 内命中的 identity 范围
    if "pident" in sub.columns and len(sub) > 0:
        sim_vmin = float(sub["pident"].min())
        sim_vmax = float(sub["pident"].max())
    else:
        sim_vmin, sim_vmax = 0.0, 100.0

    gv = GenomeViz(
        fig_width=PLOT_FIG_WIDTH,
        fig_track_height=PLOT_FIG_TRACK_HEIGHT,
        track_align_type=PLOT_TRACK_ALIGN_TYPE,
        feature_track_ratio=PLOT_FEATURE_TRACK_RATIO,
        link_track_ratio=PLOT_LINK_TRACK_RATIO,
        theme=PLOT_THEME,
        show_axis=PLOT_SHOW_AXIS,
    )

    # -------------------------
    # SprB 参考轨道
    # -------------------------
    sprb_track = gv.add_feature_track(
        "Fj SprB",
        segments=SPRB_LENGTH,
        labelsize=TRACK_LABEL_SIZE,
        labelmargin=TRACK_LABEL_MARGIN,
        align_label=TRACK_ALIGN_LABEL,
        line_kws=TRACK_LINE_KWS,
    )

    for _, row in sprb_df.sort_values("order_index").iterrows():
        fam = row["family"]
        label = short_module_name(row["module"]) if DRAW_LABELS else ""

        sprb_track.add_feature(
            int(row["start"]),
            int(row["end"]),
            label=label,
            plotstyle=FEATURE_PLOTSTYLE,
            fc=color_map.get(fam, "#cccccc"),
            ec="black",
            lw=FEATURE_LINEWIDTH,
            text_kws=dict(size=FEATURE_LABELSIZE, rotation=FEATURE_TEXT_ROTATION),
        )

    # -------------------------
    # target 轨道
    # -------------------------
    target_track = gv.add_feature_track(
        target_id,
        segments=target_len,
        labelsize=TRACK_LABEL_SIZE,
        labelmargin=TRACK_LABEL_MARGIN,
        align_label=TRACK_ALIGN_LABEL,
        line_kws=TRACK_LINE_KWS,
    )

    for _, row in sub.iterrows():
        fam = row["query_family"]
        label = fam if (DRAW_LABELS and SHOW_TARGET_LABELS) else ""

        target_track.add_feature(
            int(row["tmin"]),
            int(row["tmax"]),
            label=label,
            plotstyle=FEATURE_PLOTSTYLE,
            fc=color_map.get(fam, "#cccccc"),
            ec="black",
            lw=FEATURE_LINEWIDTH,
            text_kws=dict(size=FEATURE_LABELSIZE, rotation=FEATURE_TEXT_ROTATION),
        )

    # -------------------------
    # 连线
    # -------------------------
    for _, hit in sub.iterrows():
        fam = hit["query_family"]
        ref_rows = build_ref_rows_for_hit(sprb_df, hit)

        hit_alpha = map_similarity_to_alpha(
            float(hit["pident"]) if "pident" in hit else None,
            sim_vmin,
            sim_vmax,
        )

        for _, ref in ref_rows.iterrows():
            gv.add_link(
                ("Fj SprB", int(ref["start"]), int(ref["end"])),
                (target_id, int(hit["tmin"]), int(hit["tmax"])),
                color=color_map.get(fam, "#999999"),
                alpha=hit_alpha,
                size=LINK_SIZE,
                curve=LINK_CURVE,
            )

    gv.set_scale_bar()

    safe_target = re.sub(r"[^A-Za-z0-9_.-]+", "_", target_id)

    if PLOT_SAVE_PNG:
        gv.savefig(outdir / f"{output_prefix}_{safe_target}.png", dpi=PLOT_DPI)
    if PLOT_SAVE_PDF:
        gv.savefig(outdir / f"{output_prefix}_{safe_target}.pdf", dpi=PLOT_DPI)

    if PLOT_SAVE_HTML:
        try:
            fig = gv.plotfig(dpi=PLOT_DPI, fast_render=False)
            gv.savefig_html(
                outdir / f"{output_prefix}_{safe_target}.html",
                figure=fig,
            )
        except Exception as e:
            log(f"[WARN] html export failed for {target_id}: {e}")

def run_plotting(candidate_df: pd.DataFrame, hits_df: pd.DataFrame, output_dir: Path) -> None:
    sprb_df = load_sprb_module_table_compatible(
        module_table_path=Path(SPRB_MODULE_TABLE),
        cluster_assignments_path=Path(CLUSTER_ASSIGNMENTS_TSV),
        cluster_summary_path=Path(CLUSTER_SUMMARY_TSV),
    )

    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    if PLOT_TARGET_MODE == "top_candidates":
        if "pass_candidate_filter" in candidate_df.columns:
            plot_targets = candidate_df[candidate_df["pass_candidate_filter"] == True]["target"].head(TOP_N_CANDIDATES).tolist()
        elif "clustered" in candidate_df.columns:
            plot_targets = candidate_df[candidate_df["clustered"] == True]["target"].head(TOP_N_CANDIDATES).tolist()
        else:
            plot_targets = candidate_df["target"].head(TOP_N_CANDIDATES).tolist()
        if len(plot_targets) < TOP_N_CANDIDATES:
            extra = candidate_df["target"].head(TOP_N_CANDIDATES).tolist()
            merged = []
            for x in plot_targets + extra:
                if x not in merged:
                    merged.append(x)
            plot_targets = merged[:TOP_N_CANDIDATES]
    elif PLOT_TARGET_MODE == "specified":
        plot_targets = SPECIFIED_TARGETS
    else:
        raise ValueError(f"Unsupported PLOT_TARGET_MODE: {PLOT_TARGET_MODE}")

    for target_id in plot_targets:
        log(f"[INFO] plotting architecture (pyGenomeViz): {target_id}")
        plot_candidate_architecture_pygenomeviz(
            sprb_df=sprb_df,
            hits_df=hits_df,
            target_id=target_id,
            outdir=plot_dir
        )

# =========================================================
# 主流程
# =========================================================

def main():
    query_fasta = Path(QUERY_FASTA)
    target_fasta = Path(TARGET_FASTA)
    output_dir = Path(OUTPUT_DIR)
    logs_dir = output_dir / "logs"
    tmp_dir = output_dir / "tmp_mmseqs"

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    if not query_fasta.exists():
        raise FileNotFoundError(f"Query fasta not found: {query_fasta}")
    if not target_fasta.exists():
        raise FileNotFoundError(f"Target fasta not found: {target_fasta}")

    query_db = output_dir / "queryDB"
    target_db = output_dir / "targetDB"
    result_db = output_dir / "resultDB"

    raw_out = output_dir / "mmseqs_raw.tsv"
    filtered_out = output_dir / "mmseqs_filtered.tsv"
    plot_filtered_out = output_dir / "mmseqs_plot_filtered.tsv"
    summary_out = output_dir / "protein_hit_summary.tsv"
    summary_filtered_out = output_dir / "protein_hit_summary_filtered.tsv"
    candidate_out = output_dir / "sprb_like_candidates.tsv"
    family_summary_out = output_dir / "query_family_summary.tsv"
    by_target_dir = output_dir / "by_target"

    prepare_for_rerun(
        query_db=query_db,
        target_db=target_db,
        result_db=result_db,
        tmp_dir=tmp_dir,
        raw_out=raw_out,
        filtered_out=filtered_out,
        plot_filtered_out=plot_filtered_out,
        summary_out=summary_out,
        summary_filtered_out=summary_filtered_out,
        candidate_out=candidate_out,
        family_summary_out=family_summary_out,
        by_target_dir=by_target_dir,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    log("[INFO] building / checking MMseqs query database")
    ensure_mmseqs_db(query_fasta, query_db)

    log("[INFO] building / checking MMseqs target database")
    ensure_mmseqs_db(target_fasta, target_db)

    log("[INFO] running mmseqs search")
    run_mmseqs_search(query_db, target_db, result_db, tmp_dir)

    log("[INFO] exporting mmseqs results")
    export_mmseqs_tsv(query_db, target_db, result_db, raw_out)

    log("[INFO] parsing mmseqs results")
    raw_df = load_mmseqs_table(raw_out)
    filtered_df = filter_hits(raw_df)
    plot_df = build_plot_df(filtered_df)
    summary_df = summarize_by_protein(filtered_df)
    summary_with_filter_df = filter_candidate_proteins(summary_df)
    candidate_df = summary_with_filter_df[summary_with_filter_df["pass_candidate_filter"]].copy()
    family_summary_df = summarize_by_family(filtered_df)
    candidate_targets = set(candidate_df["target"].tolist())
    by_target_hits_df = plot_df[plot_df["target"].isin(candidate_targets)].copy()

    filtered_df.to_csv(filtered_out, sep="\t", index=False)
    plot_df.to_csv(plot_filtered_out, sep="\t", index=False)
    summary_df.to_csv(summary_out, sep="\t", index=False)
    summary_with_filter_df.to_csv(summary_filtered_out, sep="\t", index=False)
    candidate_df.to_csv(candidate_out, sep="\t", index=False)
    family_summary_df.to_csv(family_summary_out, sep="\t", index=False)
    write_target_reports(by_target_hits_df, candidate_df, output_dir)

    write_run_info(logs_dir / "run_info.txt")

    if ENABLE_PLOTTING:
        try:
            run_plotting(candidate_df, plot_df, output_dir)
        except FileNotFoundError as e:
            print(f"[WARN] plotting skipped: {e}")

    print(f"[OK] query fasta: {query_fasta}")
    print(f"[OK] target fasta: {target_fasta}")
    print(f"[OK] output dir: {output_dir}")
    print(f"[OK] results:")
    print(f"  - mmseqs_raw.tsv")
    print(f"  - mmseqs_filtered.tsv")
    print(f"  - mmseqs_plot_filtered.tsv")
    print(f"  - protein_hit_summary.tsv")
    print(f"  - protein_hit_summary_filtered.tsv")
    print(f"  - sprb_like_candidates.tsv")
    print(f"  - query_family_summary.tsv")
    print(f"  - by_target/")
    if ENABLE_PLOTTING:
        print(f"  - plots/")

if __name__ == "__main__":
    main()
