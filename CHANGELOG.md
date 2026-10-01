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

### `Changed`

- Mixed sites in the Illumina consensus sequence are now masked as `N` by default, consistent with the Nanopore workflow. Set `--allow_iupac_codes` to restore the previous IUPAC ambiguity code behaviour.
- `--min_allele_frequency`, `--min_mask_allele_frequency` and `--min_minor_allele_count` now control variant filtering in both the Nanopore and Illumina workflows, with explicit defaults of `0.8`, `0.2` and `4` respectively.
- Nanopore variants with an allele frequency between `0.6` and `0.8` are now masked as `N` rather than applied, and variants between `0.1` and `0.2` are now discarded rather than masked (previously the artic minion defaults of `0.6` / `0.1` were used).
- Illumina variants now require `4` (previously `10`) supporting reads to be called, and variants with an allele frequency of exactly `--min_allele_frequency` are now applied rather than treated as mixed, consistent with artic minion.

### `Fixed`

- Nanopore samples with an explicit FASTQ path and a `barcode` were silently dropped, explicit paths now always take precedence over `barcode` fuzzy matching.
- Fixed mismatched header / column order in the example `assets/samplesheet.csv`.

### `Dependencies`

### `Deprecated`

### `Removed`

- The Illumina-specific `--lower_ambiguity_frequency`, `--upper_ambiguity_frequency` and `--min_ambiguity_count` parameters have been removed, use `--min_mask_allele_frequency`, `--min_allele_frequency` and `--min_minor_allele_count` respectively instead.
