# artic-network/amplicon-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.12.0dev - [date]

### `Added`

- Nanopore samples can now be provided as a single FASTQ file in the `fastq_1` samplesheet column, as an alternative to a directory of FASTQ files in `fastq_directory`. Reads are still filtered by `artic guppyplex` in both cases.
- Uncompressed FASTQ files (`.fastq` / `.fq`) are now accepted in the `fastq_1` and `fastq_2` columns for both Nanopore and Illumina samples.
- A warning is now raised when an Illumina sample has a `fastq_directory` value (which is ignored) or a Nanopore sample has a `fastq_2` value (which is ignored).
- A Nanopore sample which provides both `fastq_directory` and `fastq_1` now fails samplesheet validation with an informative error.
- New `--allow_iupac_codes` parameter (default `false`) to represent mixed sites in the Illumina consensus sequence with IUPAC ambiguity codes instead of `N`.
- A Pathoplexus submission folder is now prepared per primer scheme by default (turn off with the new `--skip_pathoplexus_template` parameter).
  - Each folder contains a virus-agnostic metadata template (XLSX and TSV), a FASTA with matching IDs, and the dehosted reads with FASTQ headers stripped to the read ID.
  - Everything that can be derived from the run is filled in automatically, including sequencing date and flow cell type from Nanopore read headers, QC results and dropout amplicons.
  - Not designed for primer panels or multi-pathogen schemes, where the reference contains several sequences that are not segments of one genome (segmented genomes and reference-selecting schemes are supported).
- The per-sample QC TSV (`<SAMPLE>_amplicon-nf_qc-report.tsv`) has a new `amplicon_dropouts` column listing the amplicons below `--min_coverage_depth`.
- `align_trim` QC reports are now published per sample for both Nanopore and Illumina samples: `<SAMPLE>.summary.tsv`, `<SAMPLE>.amplicons.tsv`, `<SAMPLE>.primer_pairs.tsv`, `<SAMPLE>.split_reads.tsv` and `<SAMPLE>.overlaps.tsv`. They show how many reads were kept or removed and why, per-amplicon read fates, the primer pairs observed, split (chimeric) alignments and the concordance between overlapping amplicons.
- The run report has a new Primer Trimming QC tab. It summarises the `align_trim` QC reports for every sample in the scheme: read fates, per-amplicon fates, primer pair types, split read joins and amplicon overlaps. Values are shown without pass/fail thresholds.
- The run report can now show run-level warnings in a red box at the top. A warning is raised when pre- and post-normalisation amplicon coverage disagree by more than 1%, for example when `--normalise_depth` is lower than `--min_coverage_depth`.
- The run QC TSV (`<SCHEME_NAME>_qc_results.tsv`) has new `pre_norm_amplicon_coverage` and `post_norm_amplicon_coverage` columns, appended at the end. They give the percentage of positions inside amplicon inserts (primers excluded) with depth ≥ `--min_coverage_depth`.

### `Changed`

- Mixed sites in the Illumina consensus sequence are now masked as `N` by default, consistent with the Nanopore workflow. Set `--allow_iupac_codes` to restore the previous IUPAC ambiguity code behaviour.
- `--min_allele_frequency`, `--min_mask_allele_frequency` and `--min_minor_allele_count` now control variant filtering in both the Nanopore and Illumina workflows, with explicit defaults of `0.8`, `0.2` and `4` respectively.
- Nanopore variants with an allele frequency between `0.6` and `0.8` are now masked as `N` rather than applied, and variants between `0.1` and `0.2` are now discarded rather than masked (previously the artic minion defaults of `0.6` / `0.1` were used).
- Illumina variants now require `4` (previously `10`) supporting reads to be called, and variants with an allele frequency of exactly `--min_allele_frequency` are now applied rather than treated as mixed, consistent with artic minion.
- The run report QC table has been simplified.
  - The separate pre/post-normalisation coverage columns are replaced by a single "Raw Coverage" column (pre-normalisation, whole reference) and an "Amplicon Coverage" column (amplicon inserts only, post-normalisation when normalising).
  - The Primer Scheme column has been removed, since the scheme is already shown in the report header.
  - The run QC TSV still contains all the previous columns.

### `Fixed`

- The run and sample reports now read the `mean_fragment_depth` column written by `align_trim` 2.0, which replaces `mean_depth`. Previously both reports would fail.
- Nanopore samples with an explicit FASTQ path and a `barcode` were silently dropped, explicit paths now always take precedence over `barcode` fuzzy matching.
- Fixed mismatched header / column order in the example `assets/samplesheet.csv`.
- Illumina samples now use their own reference when a primer scheme picks a reference variant per sample (refselect schemes such as `artic-pan-ebola`). Previously the reference, BWA index and FASTA index were deduplicated per scheme, so every Illumina sample of the scheme was aligned, variant called and given a consensus against the first sample's variant. These are now deduplicated per distinct reference.
- `--min_mapping_quality` is now passed to freebayes in the Illumina workflow. Previously freebayes applied its own default (`1`), discarding MAPQ 0 reads kept by `align_trim`, which caused well-covered repetitive regions (e.g. mpox inverted terminal repeats) to be depth masked when `--min_mapping_quality 0` was set.

### `Dependencies`

| Dependency                     | Old version | New version |
| ------------------------------ | ----------- | ----------- |
| `artic`                        | 1.11.2      | 1.12.0      |
| `align_trim`                   | 1.2.0       | 2.0.0       |
| `samtools` (`ARTIC_ALIGNTRIM`) | 1.23.1      | 1.25        |

### `Deprecated`

### `Removed`

- The Illumina-specific `--lower_ambiguity_frequency`, `--upper_ambiguity_frequency` and `--min_ambiguity_count` parameters have been removed, use `--min_mask_allele_frequency`, `--min_allele_frequency` and `--min_minor_allele_count` respectively instead.
