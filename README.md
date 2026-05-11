# SprB Genome Scan

This repository now includes a modular SprB-like protein genome scan pipeline extracted from the prototype `sprb_module_mmseqs_scanner.py` while preserving its filtering logic, output fields, and plotting flow.

The workflow is now split into two stages:

- `run_sprb_scan.py`: HPC-friendly scan stage
- `run_sprb_plot.py`: local plotting stage

## Layout

- `src/config.py`: load JSON/YAML config into typed parameter objects
- `src/load_panel.py`: load standard SprB module tables or domain table plus cluster mapping
- `src/scanner.py`: run MMseqs2 and build hit-level DataFrames
- `src/hit_processing.py`: compute protein-level hit summaries
- `src/candidate_filter.py`: apply SprB-like candidate thresholds
- `src/plotting.py`: render pyGenomeViz architecture plots
- `configs/scan_config.json`: scan-stage configuration
- `configs/plot_config.json`: plot-stage configuration
- `run_sprb_scan.py`: HPC-friendly scan runner
- `run_sprb_plot.py`: local plotting runner

## Stage 1: Scan

Update paths in [configs/scan_config.json](/Users/shulei/PycharmProjects/Biopython/sprB/sprb_genome_scan/configs/scan_config.json), then run:

```bash
python3 run_sprb_scan.py
```

`run_sprb_scan.py` reads its default config directly from [configs/scan_config.json](/Users/shulei/PycharmProjects/Biopython/sprB/sprb_genome_scan/configs/scan_config.json). This stage is designed for HPC use and does not depend on plotting settings.

## Target Inputs

The scan stage supports two native input styles:

- Single FASTA: set `target_input` to one `.faa` file and use `target_mode: "single_fasta"` or `target_mode: "auto"`
- NCBI genome bundle directory: set `target_input` to a directory like `ncbi_gliding_genomes/` and use `target_mode: "ncbi_genome_dir"`

For NCBI-style bundle directories, the scanner will automatically walk each genome subdirectory and pick the first file matching `target_glob` such as `*_protein.faa` or `*_protein.faa.gz`.

## Stage 2: Plot

After copying the scan outputs back to your local machine, update [configs/plot_config.json](/Users/shulei/PycharmProjects/Biopython/sprB/sprb_genome_scan/configs/plot_config.json), then run:

```bash
python3 run_sprb_plot.py
```

`run_sprb_plot.py` only needs the scan result tables plus the SprB module reference files. It does not touch MMseqs or the original genome FASTA inputs.

## Outputs

The scan stage saves the compact result files needed for downstream work:

- `hits.tsv`
- `protein_hit_summary.tsv`
- `sprb_like_candidates.tsv`
- per-genome subdirectories with their own `hits.tsv`, `protein_hit_summary.tsv`, and `sprb_like_candidates.tsv`

The plot stage writes:

- `plots/` under each plotted output directory
- `logs/plot_run_info.txt`

The scan stage writes `logs/scan_run_info.txt`, and the plot stage writes `logs/plot_run_info.txt`.
