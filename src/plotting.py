from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import PlotStyleConfig
from .load_panel import load_sprb_module_table
from .utils import safe_filename, short_module_name


def build_family_color_map(families: list[str]) -> dict[str, str]:
    unique_fams = sorted(set(families))
    base_colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1",
        "#9c755f", "#ff9da7", "#bab0ab", "#86bc86", "#499894",
    ]
    return {family: base_colors[i % len(base_colors)] for i, family in enumerate(unique_fams)}


def build_ref_rows_for_block(sprb_df: pd.DataFrame, block_row: pd.Series, plot_cfg: PlotStyleConfig) -> pd.DataFrame:
    family = block_row["best_family"]
    member = block_row["best_member"]
    if plot_cfg.link_mode == "all":
        return sprb_df[sprb_df["family"] == family]
    if plot_cfg.link_mode == "best":
        ref_rows = sprb_df[
            (sprb_df["family"].astype(str) == str(family)) &
            (sprb_df["module"].astype(str) == str(member))
        ]
        if ref_rows.empty:
            member_domain = f"domain_{member}" if not str(member).startswith("domain_") else str(member)
            ref_rows = sprb_df[
                (sprb_df["family"].astype(str) == str(family)) &
                (sprb_df["module"].astype(str) == member_domain)
            ]
        if ref_rows.empty:
            ref_rows = sprb_df[sprb_df["family"] == family].head(1)
        return ref_rows
    raise ValueError(f"Unsupported PLOT_LINK_MODE: {plot_cfg.link_mode}")


def map_similarity_to_alpha(
    value: float,
    vmin: float,
    vmax: float,
    alpha_min: float,
    alpha_max: float,
) -> float:
    if pd.isna(value):
        return alpha_min
    if vmax <= vmin:
        return alpha_max
    scaled = (value - vmin) / (vmax - vmin)
    scaled = max(0.0, min(1.0, scaled))
    return alpha_min + scaled * (alpha_max - alpha_min)


def plot_candidate_architecture(
    sprb_df: pd.DataFrame,
    blocks_df: pd.DataFrame,
    target_id: str,
    outdir: str | Path,
    plot_cfg: PlotStyleConfig,
    output_prefix: str = "sprb_candidate_architecture",
) -> None:
    from pygenomeviz import GenomeViz

    sub = blocks_df[blocks_df["target"] == target_id].copy()
    if sub.empty:
        return

    sub = sub.sort_values("block_start").copy()
    target_len = int(sub["tlen"].max())
    color_map = build_family_color_map(list(sprb_df["family"]) + list(sub["best_family"]))
    sim_vmin = float(sub["mean_pident"].min()) if "mean_pident" in sub.columns and len(sub) > 0 else 0.0
    sim_vmax = float(sub["mean_pident"].max()) if "mean_pident" in sub.columns and len(sub) > 0 else 100.0

    gv = GenomeViz(
        fig_width=plot_cfg.fig_width,
        fig_track_height=plot_cfg.fig_track_height,
        track_align_type=plot_cfg.track_align_type,
        feature_track_ratio=plot_cfg.feature_track_ratio,
        link_track_ratio=plot_cfg.link_track_ratio,
        theme=plot_cfg.theme,
        show_axis=plot_cfg.show_axis,
    )
    sprb_track = gv.add_feature_track(
        "Fj SprB",
        segments=plot_cfg.sprb_length,
        labelsize=plot_cfg.track_label_size,
        labelmargin=plot_cfg.track_label_margin,
        align_label=plot_cfg.track_align_label,
        line_kws=plot_cfg.track_line_kws,
    )
    for _, row in sprb_df.sort_values("order_index").iterrows():
        sprb_track.add_feature(
            int(row["start"]),
            int(row["end"]),
            label=short_module_name(row["module"]) if plot_cfg.draw_labels else "",
            plotstyle=plot_cfg.feature_plotstyle,
            fc=color_map.get(row["family"], "#cccccc"),
            ec="black",
            lw=plot_cfg.feature_linewidth,
            text_kws=dict(size=plot_cfg.feature_labelsize, rotation=plot_cfg.feature_text_rotation),
        )

    target_track = gv.add_feature_track(
        target_id,
        segments=target_len,
        labelsize=plot_cfg.track_label_size,
        labelmargin=plot_cfg.track_label_margin,
        align_label=plot_cfg.track_align_label,
        line_kws=plot_cfg.track_line_kws,
    )
    for _, row in sub.iterrows():
        family_label = str(row["best_family"])
        target_track.add_feature(
            int(row["block_start"]),
            int(row["block_end"]),
            label=family_label if (plot_cfg.draw_labels and plot_cfg.show_target_labels) else "",
            plotstyle=plot_cfg.feature_plotstyle,
            fc=color_map.get(family_label, "#cccccc"),
            ec="black",
            lw=plot_cfg.feature_linewidth,
            text_kws=dict(size=plot_cfg.feature_labelsize, rotation=plot_cfg.feature_text_rotation),
        )

    for _, block in sub.iterrows():
        family_label = str(block["best_family"])
        alpha = map_similarity_to_alpha(
            float(block["mean_pident"]) if "mean_pident" in block else None,
            sim_vmin,
            sim_vmax,
            plot_cfg.link_alpha_min,
            plot_cfg.link_alpha_max,
        )
        ref_rows = build_ref_rows_for_block(sprb_df, block, plot_cfg)
        for _, ref in ref_rows.iterrows():
            gv.add_link(
                ("Fj SprB", int(ref["start"]), int(ref["end"])),
                (target_id, int(block["block_start"]), int(block["block_end"])),
                color=color_map.get(family_label, "#999999"),
                alpha=alpha,
                size=plot_cfg.link_size,
                curve=plot_cfg.link_curve,
            )

    gv.set_scale_bar()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    safe_target = safe_filename(target_id)
    if plot_cfg.save_png:
        gv.savefig(outdir / f"{output_prefix}_{safe_target}.png", dpi=plot_cfg.dpi)
    if plot_cfg.save_pdf:
        gv.savefig(outdir / f"{output_prefix}_{safe_target}.pdf", dpi=plot_cfg.dpi)
    if plot_cfg.save_html:
        fig = gv.plotfig(dpi=plot_cfg.dpi, fast_render=False)
        gv.savefig_html(outdir / f"{output_prefix}_{safe_target}.html", figure=fig)


def _select_plot_targets(plot_source_df: pd.DataFrame, plot_cfg: PlotStyleConfig) -> list[str]:
    if plot_cfg.target_mode == "top_candidates":
        if "pass_candidate_filter" in plot_source_df.columns and plot_source_df["pass_candidate_filter"].any():
            plot_targets = plot_source_df[plot_source_df["pass_candidate_filter"] == True]["target"].head(plot_cfg.top_n_candidates).tolist()
        else:
            plot_targets = plot_source_df["target"].head(plot_cfg.top_n_candidates).tolist()
        if len(plot_targets) < plot_cfg.top_n_candidates:
            extra = plot_source_df["target"].head(plot_cfg.top_n_candidates).tolist()
            merged = []
            for target in plot_targets + extra:
                if target not in merged:
                    merged.append(target)
            plot_targets = merged[:plot_cfg.top_n_candidates]
        return plot_targets
    if plot_cfg.target_mode == "specified":
        return plot_cfg.specified_targets
    raise ValueError(f"Unsupported PLOT_TARGET_MODE: {plot_cfg.target_mode}")


def run_plotting(plot_source_df: pd.DataFrame, blocks_df: pd.DataFrame, plot_cfg: PlotStyleConfig, output_dir: str | Path) -> list[str]:
    sprb_df = load_sprb_module_table(
        module_table_path=plot_cfg.sprb_module_table,
        cluster_assignments_path=plot_cfg.cluster_assignments_tsv or None,
        cluster_summary_path=plot_cfg.cluster_summary_tsv or None,
    )
    plot_targets = _select_plot_targets(plot_source_df, plot_cfg)
    plot_dir = Path(output_dir) / "plots"
    for target_id in plot_targets:
        plot_candidate_architecture(sprb_df, blocks_df, target_id, plot_dir, plot_cfg)
    return plot_targets


def run_plotting_for_result_dir(result_dir: str | Path, output_root: str | Path, plot_cfg: PlotStyleConfig) -> list[str]:
    result_dir = Path(result_dir)
    blocks_df = pd.read_csv(result_dir / "blocks.tsv", sep="\t")
    summary_df = pd.read_csv(result_dir / "protein_hit_summary.tsv", sep="\t")
    candidate_df = pd.read_csv(result_dir / "sprb_like_candidates.tsv", sep="\t")
    plot_source_df = candidate_df if len(candidate_df) > 0 else summary_df
    destination = Path(output_root) / result_dir.name if Path(output_root) != result_dir else result_dir
    return run_plotting(plot_source_df, blocks_df, plot_cfg, destination)
