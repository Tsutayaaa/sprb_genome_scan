# SprB Genome Scan

This repository now includes a modular SprB-like protein genome scan pipeline extracted from the prototype `sprb_module_mmseqs_scanner.py` while preserving its filtering logic, output fields, and plotting flow.

## Layout

- `src/config.py`: load JSON/YAML config into typed parameter objects
- `src/load_panel.py`: load standard SprB module tables or domain table plus cluster mapping
- `src/scanner.py`: run MMseqs2 and build hit-level DataFrames
- `src/hit_processing.py`: compute protein-level hit summaries
- `src/candidate_filter.py`: apply SprB-like candidate thresholds
- `src/plotting.py`: render pyGenomeViz architecture plots
- `src/pipeline.py`: orchestrate the full run and save final outputs
- `configs/default_config.yaml`: example configuration
- `run_sprb_scan.py`: entrypoint script

## Usage

Update paths in [configs/default_config.yaml](/Users/shulei/PycharmProjects/Biopython/sprB/sprb_genome_scan/configs/default_config.yaml), then run:

```bash
python3 run_sprb_scan.py configs/default_config.yaml
```

## Outputs

Only plotting-level filtered hits and downstream summaries are saved:

- `hits.tsv`
- `protein_hit_summary.tsv`
- `sprb_like_candidates.tsv`
- `plots/`

`logs/run_info.txt` captures the effective run configuration.
