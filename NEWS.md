# Cysinger 1.0.0

*Released: 2026-02-22*

## Initial Release

### Pipeline Functions (replacing 17 Perl scripts)
- `filter_fastq_reads()` — FASTQ read filtering by MID tags and primers
- `demultiplex_reads()` — Bin assignment by MID × primer combinations
- `quality_filter_reads()` — Phred Q20 quality filtering and length trimming
- `reads_to_fasta()` — FASTQ to FASTA conversion with reverse complement
- `align_to_reference()` — Pairwise alignment via Biostrings (no external EMBOSS dependency)
- `combine_paired_reads()` — Forward/reverse read merging with overlap handling
- `build_substitution_table()` — WT vs mutant codon comparison
- `classify_mutations()` — WT/single/double/triple/multi classification
- `calculate_cys_frequency()` — Codon-level cysteine frequency counting
- `aggregate_cys_residue_level()` — Residue-level frequency aggregation
- `normalize_cys_frequency()` — Enrichment ratio normalization (sorted ÷ input)
- `run_pipeline()` — Full end-to-end pipeline wrapper

### Visualization Functions (7 publication-quality figures)
- `plot_cys_heatmap()` — Inferno-colored frequency heatmap
- `plot_cys_enrichment()` — Grouped bar plot of read counts per bin
- `plot_mutation_landscape()` — Lollipop plot with auto-labeled top hits
- `plot_epitope_profile()` — Loess-smoothed enrichment profile with region highlights
- `plot_surface_accessibility()` — Polar radar chart of accessible vs buried residues
- `plot_qc_summary()` — 4-panel QC dashboard (patchwork composite)
- `plot_mutation_class_donut()` — Donut chart of mutation classification proportions

### Utilities
- `get_codon_table()` — Standard genetic code (64 codons)
- `translate_sequence()` — Nucleotide to protein translation
- `reverse_complement_seq()` — DNA reverse complement
- `ccdB_wt_nucleotide()` / `ccdB_wt_protein()` — CcdB reference sequences
- `load_example_data()` — Built-in CcdB dataset loader
- `generate_example_normalized()` — Pre-normalized demo dataset

### Data
- Included CcdB cysteine scanning data from Najar et al. (*Structure*, 2017)
- Wild-type CcdB nucleotide sequence (306 nt) and protein sequence (101 aa)
