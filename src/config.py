from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import json


@dataclass
class MMseqsConfig:
    sensitivity: float = 7.5
    min_seq_id: float = 0.15
    evalue: float = 1e-3
    max_seqs: int = 5000
    cov_mode: int = 0
    min_coverage: float = 0.0
    alt_ali: int = 20


@dataclass
class FilterHitsConfig:
    min_bits: float = 30.0
    min_aln_len: int = 50
    min_pident: float = 15.0


@dataclass
class FilterProteinConfig:
    min_length: int = 500
    min_blocks: int = 2
    min_unique_families: int = 2
    min_unique_members: int = 2
    min_cluster_span_coverage: float = 0.20
    min_union_covered_len: int = 200
    min_union_coverage: float = 0.08
    min_linearity_score: float = 0.70
    max_order_inversions: int = 1
    require_repeat_pattern: bool = True
    repeat_min_motif_length: int = 2
    repeat_max_motif_length: int = 4
    repeat_min_occurrences: int = 2
    repeat_max_gap_blocks: int = 1
    min_repeat_count: int = 2


@dataclass
class PlotHitFilterConfig:
    min_bits: float = 35.0
    min_aln_len: int = 60
    min_pident: float = 15.0


@dataclass
class BlockMergeConfig:
    overlap_rate_threshold: float = 0.90


@dataclass
class AssemblyMetadataConfig:
    enabled: bool = True
    report_path: str = ""
    report_path_relative_to: str = "target_input"
    display_name_source: str = "organism_name"
    output_dir_source: str = "organism_plus_accession"


@dataclass
class ReferencePanelConfig:
    sprb_module_table: str = ""
    cluster_assignments_tsv: str = ""
    cluster_summary_tsv: str = ""


@dataclass
class PlotStyleConfig:
    sprb_module_table: str = ""
    cluster_assignments_tsv: str = ""
    cluster_summary_tsv: str = ""
    sprb_length: int = 6497
    top_n_candidates: int = 10
    target_mode: str = "top_candidates"
    specified_targets: list[str] = field(default_factory=list)
    link_mode: str = "all"
    save_png: bool = True
    save_pdf: bool = True
    save_html: bool = False
    fig_width: int = 18
    fig_track_height: float = 1.2
    track_align_type: str = "center"
    feature_track_ratio: float = 0.32
    link_track_ratio: float = 0.90
    theme: str = "light"
    show_axis: bool = False
    track_label_size: int = 16
    track_label_margin: float = 0.02
    track_align_label: bool = True
    track_line_kws: dict = field(default_factory=lambda: {"color": "black", "lw": 0.8})
    draw_labels: bool = True
    feature_plotstyle: str = "box"
    feature_linewidth: float = 0.6
    feature_labelsize: int = 8
    feature_text_rotation: int = 90
    show_target_labels: bool = False
    link_curve: bool = True
    link_size: float = 0.9
    link_alpha_min: float = 0.12
    link_alpha_max: float = 0.80
    dpi: int = 300


@dataclass
class ScanConfig:
    query_fasta: str
    output_dir: str
    target_fasta: str = ""
    target_input: str = ""
    target_mode: str = "auto"
    target_glob: str = "*_protein.faa"
    mmseqs_bin: str = "mmseqs"
    threads: int = 8
    force_rerun: bool = True
    verbose: bool = True
    mmseqs: MMseqsConfig = field(default_factory=MMseqsConfig)
    filter_hits: FilterHitsConfig = field(default_factory=FilterHitsConfig)
    filter_protein: FilterProteinConfig = field(default_factory=FilterProteinConfig)
    plot_filter_hits: PlotHitFilterConfig = field(default_factory=PlotHitFilterConfig)
    block_merge: BlockMergeConfig = field(default_factory=BlockMergeConfig)
    assembly_metadata: AssemblyMetadataConfig = field(default_factory=AssemblyMetadataConfig)
    reference_panel: ReferencePanelConfig = field(default_factory=ReferencePanelConfig)


@dataclass
class PlotConfig:
    results_root: str
    output_dir: str = ""
    results_mode: str = "auto"
    verbose: bool = True
    plot: PlotStyleConfig = field(default_factory=PlotStyleConfig)


def _load_yaml_or_json(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        return json.loads(text)

    try:
        import yaml  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "YAML config requested but PyYAML is not installed. Use JSON or install PyYAML."
        ) from exc
    return yaml.safe_load(text)


def load_scan_config(path: str | Path) -> ScanConfig:
    raw = _load_yaml_or_json(Path(path))
    config = ScanConfig(
        mmseqs=MMseqsConfig(**raw.get("mmseqs", {})),
        filter_hits=FilterHitsConfig(**raw.get("filter_hits", {})),
        filter_protein=FilterProteinConfig(**raw.get("filter_protein", {})),
        plot_filter_hits=PlotHitFilterConfig(**raw.get("plot_filter_hits", raw.get("plot_hits_filter", {}))),
        block_merge=BlockMergeConfig(**raw.get("block_merge", {})),
        assembly_metadata=AssemblyMetadataConfig(**raw.get("assembly_metadata", {})),
        reference_panel=ReferencePanelConfig(**raw.get("reference_panel", {})),
        **{
            key: value
            for key, value in raw.items()
            if key not in {"mmseqs", "filter_hits", "filter_protein", "plot_filter_hits", "plot_hits_filter", "plot", "assembly_metadata", "block_merge", "reference_panel"}
        },
    )
    if not config.target_input:
        config.target_input = config.target_fasta
    if not config.target_fasta:
        config.target_fasta = config.target_input
    return config


def load_plot_config(path: str | Path) -> PlotConfig:
    raw = _load_yaml_or_json(Path(path))
    config = PlotConfig(
        plot=PlotStyleConfig(**raw.get("plot", {})),
        **{key: value for key, value in raw.items() if key != "plot"},
    )
    if not config.output_dir:
        config.output_dir = config.results_root
    return config


def config_to_dict(config: object) -> dict:
    return asdict(config)
