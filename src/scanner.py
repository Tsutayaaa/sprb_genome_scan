from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import logging
import gzip

import pandas as pd

from .config import FilterHitsConfig, MMseqsConfig, PlotHitFilterConfig, ScanConfig
from .utils import extract_family, extract_member, remove_prefix_artifacts, run_command


RAW_COLUMNS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen",
]


@dataclass
class ScanPaths:
    query_db: Path
    target_db: Path
    result_db: Path
    tmp_dir: Path
    raw_export_tsv: Path
    hits_out: Path
    summary_out: Path
    candidates_out: Path
    logs_dir: Path
    plots_dir: Path


@dataclass
class TargetSpec:
    genome_name: str
    fasta_path: Path
    output_dir: Path


def build_scan_paths(output_dir: str | Path) -> ScanPaths:
    output_dir = Path(output_dir)
    tmp_dir = output_dir / "tmp_mmseqs"
    return ScanPaths(
        query_db=tmp_dir / "queryDB",
        target_db=tmp_dir / "targetDB",
        result_db=tmp_dir / "resultDB",
        tmp_dir=tmp_dir,
        raw_export_tsv=tmp_dir / "mmseqs_raw.tsv",
        hits_out=output_dir / "hits.tsv",
        summary_out=output_dir / "protein_hit_summary.tsv",
        candidates_out=output_dir / "sprb_like_candidates.tsv",
        logs_dir=output_dir / "logs",
        plots_dir=output_dir / "plots",
    )


def _strip_known_suffixes(path: Path) -> str:
    name = path.name
    for suffix in [".faa.gz", ".fasta.gz", ".fa.gz", ".faa", ".fasta", ".fa"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def discover_targets(config: ScanConfig) -> list[TargetSpec]:
    target_input = Path(config.target_input or config.target_fasta)
    if not target_input.exists():
        raise FileNotFoundError(f"Target input not found: {target_input}")

    mode = config.target_mode.lower()
    if mode not in {"auto", "single_fasta", "ncbi_genome_dir"}:
        raise ValueError(f"Unsupported target_mode: {config.target_mode}")

    output_root = Path(config.output_dir)
    if target_input.is_file():
        genome_name = _strip_known_suffixes(target_input)
        return [TargetSpec(genome_name=genome_name, fasta_path=target_input, output_dir=output_root)]

    fasta_paths: list[Path] = []
    if mode in {"auto", "ncbi_genome_dir"}:
        for child in sorted(target_input.iterdir()):
            if not child.is_dir():
                continue
            matches = sorted(child.glob(config.target_glob))
            if not matches:
                matches = sorted(child.glob(config.target_glob + ".gz"))
            if matches:
                fasta_paths.append(matches[0])

    if not fasta_paths and mode == "auto":
        fasta_paths = sorted(target_input.glob(config.target_glob))
        if not fasta_paths:
            fasta_paths = sorted(target_input.glob(config.target_glob + ".gz"))

    if not fasta_paths:
        raise FileNotFoundError(
            f"No target FASTA files matching '{config.target_glob}' were found under: {target_input}"
        )

    targets = []
    for fasta_path in fasta_paths:
        genome_name = fasta_path.parent.name if fasta_path.parent != target_input else _strip_known_suffixes(fasta_path)
        targets.append(
            TargetSpec(
                genome_name=genome_name,
                fasta_path=fasta_path,
                output_dir=output_root / genome_name,
            )
        )
    return targets


def _materialize_target_fasta(target_path: Path, tmp_dir: Path) -> Path:
    if target_path.suffix != ".gz":
        return target_path
    materialized = tmp_dir / _strip_known_suffixes(target_path)
    with gzip.open(target_path, "rt", encoding="utf-8") as src, materialized.open("w", encoding="utf-8") as dst:
        dst.write(src.read())
    return materialized


def prepare_for_rerun(paths: ScanPaths, force_rerun: bool, logger: logging.Logger | None = None) -> None:
    if not force_rerun:
        return
    if logger:
        logger.info("FORCE_RERUN=True, removing old MMseqs artifacts")
    remove_prefix_artifacts(paths.query_db)
    remove_prefix_artifacts(paths.target_db)
    remove_prefix_artifacts(paths.result_db)
    if paths.tmp_dir.exists():
        shutil.rmtree(paths.tmp_dir)
    for file_path in [paths.hits_out, paths.summary_out, paths.candidates_out]:
        if file_path.exists():
            file_path.unlink()
    if paths.plots_dir.exists():
        shutil.rmtree(paths.plots_dir)


def ensure_mmseqs_db(fasta: Path, db_path: Path, mmseqs_bin: str, logger: logging.Logger | None = None, force_rerun: bool = False) -> None:
    if db_path.with_suffix(".dbtype").exists() and not force_rerun:
        if logger:
            logger.info("MMseqs db already exists: %s", db_path)
        return
    run_command([mmseqs_bin, "createdb", str(fasta), str(db_path)], logger=logger)


def run_mmseqs_search(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    tmp_dir: Path,
    mmseqs_bin: str,
    threads: int,
    mmseqs_cfg: MMseqsConfig,
    logger: logging.Logger | None = None,
    force_rerun: bool = False,
) -> None:
    if result_db.with_suffix(".dbtype").exists() and not force_rerun:
        if logger:
            logger.info("MMseqs result db already exists: %s", result_db)
        return

    cmd = [
        mmseqs_bin, "search",
        str(query_db), str(target_db), str(result_db), str(tmp_dir),
        "--threads", str(threads),
        "-s", str(mmseqs_cfg.sensitivity),
        "-e", str(mmseqs_cfg.evalue),
        "--min-seq-id", str(mmseqs_cfg.min_seq_id),
        "--max-seqs", str(mmseqs_cfg.max_seqs),
        "--cov-mode", str(mmseqs_cfg.cov_mode),
        "-c", str(mmseqs_cfg.min_coverage),
        "--alt-ali", str(mmseqs_cfg.alt_ali),
    ]
    run_command(cmd, logger=logger)


def export_mmseqs_tsv(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    out_tsv: Path,
    mmseqs_bin: str,
    logger: logging.Logger | None = None,
    force_rerun: bool = False,
) -> None:
    if out_tsv.exists() and not force_rerun:
        if logger:
            logger.info("MMseqs TSV already exists: %s", out_tsv)
        return

    outfmt = ",".join(RAW_COLUMNS)
    cmd = [
        mmseqs_bin, "convertalis",
        str(query_db), str(target_db), str(result_db), str(out_tsv),
        "--format-output", outfmt,
    ]
    run_command(cmd, logger=logger)


def load_mmseqs_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=RAW_COLUMNS)

    df = pd.read_csv(path, sep="\t", header=None, names=RAW_COLUMNS, dtype=str)
    if len(df) > 0 and str(df.iloc[0]["query"]) == "query":
        df = df.iloc[1:].copy()

    numeric_cols = [
        "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
        "tstart", "tend", "evalue", "bits", "qlen", "tlen",
    ]
    for column in numeric_cols:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=["fident", "alnlen", "bits", "tstart", "tend", "tlen"]).copy()
    df["pident"] = df["fident"] * 100.0
    df["tmin"] = df[["tstart", "tend"]].min(axis=1)
    df["tmax"] = df[["tstart", "tend"]].max(axis=1)
    df["strand"] = df.apply(lambda row: "+" if row["tend"] >= row["tstart"] else "-", axis=1)
    df["query_family"] = df["query"].map(extract_family)
    df["query_member"] = df["query"].map(extract_member)
    return df


def filter_hits(df: pd.DataFrame, filter_cfg: FilterHitsConfig) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return df[
        (df["pident"] >= filter_cfg.min_pident) &
        (df["alnlen"] >= filter_cfg.min_aln_len) &
        (df["bits"] >= filter_cfg.min_bits)
    ].copy()


def build_plot_df(df: pd.DataFrame, filter_cfg: PlotHitFilterConfig) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return df[
        (df["bits"] >= filter_cfg.min_bits) &
        (df["alnlen"] >= filter_cfg.min_aln_len) &
        (df["pident"] >= filter_cfg.min_pident)
    ].copy()


def run_scan_for_target(
    config: ScanConfig,
    target_spec: TargetSpec,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    query_fasta = Path(config.query_fasta)
    if not query_fasta.exists():
        raise FileNotFoundError(f"Query fasta not found: {query_fasta}")
    if not target_spec.fasta_path.exists():
        raise FileNotFoundError(f"Target fasta not found: {target_spec.fasta_path}")

    paths = build_scan_paths(target_spec.output_dir)
    for path in [target_spec.output_dir, paths.logs_dir, paths.tmp_dir]:
        path.mkdir(parents=True, exist_ok=True)

    prepare_for_rerun(paths, config.force_rerun, logger=logger)
    for path in [target_spec.output_dir, paths.logs_dir, paths.tmp_dir]:
        path.mkdir(parents=True, exist_ok=True)

    target_fasta = _materialize_target_fasta(target_spec.fasta_path, paths.tmp_dir)
    ensure_mmseqs_db(query_fasta, paths.query_db, config.mmseqs_bin, logger=logger, force_rerun=config.force_rerun)
    ensure_mmseqs_db(target_fasta, paths.target_db, config.mmseqs_bin, logger=logger, force_rerun=config.force_rerun)
    run_mmseqs_search(
        paths.query_db,
        paths.target_db,
        paths.result_db,
        paths.tmp_dir,
        config.mmseqs_bin,
        config.threads,
        config.mmseqs,
        logger=logger,
        force_rerun=config.force_rerun,
    )
    export_mmseqs_tsv(
        paths.query_db,
        paths.target_db,
        paths.result_db,
        paths.raw_export_tsv,
        config.mmseqs_bin,
        logger=logger,
        force_rerun=config.force_rerun,
    )
    raw_df = load_mmseqs_table(paths.raw_export_tsv)
    if paths.raw_export_tsv.exists():
        paths.raw_export_tsv.unlink()
    filtered_df = filter_hits(raw_df, config.filter_hits)
    plot_df = build_plot_df(filtered_df, config.plot_filter_hits)
    for df in (filtered_df, plot_df):
        if df.empty:
            df["genome_name"] = pd.Series(dtype=str)
            df["target_fasta"] = pd.Series(dtype=str)
        else:
            df["genome_name"] = target_spec.genome_name
            df["target_fasta"] = str(target_spec.fasta_path)
    return filtered_df, plot_df
