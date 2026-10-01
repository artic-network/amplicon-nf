# artic-network/amplicon-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.12.0dev - [date]

### `Added`

- Nanopore samples can now be provided as a single FASTQ file in the `fastq_1` samplesheet column, as an alternative to a directory of FASTQ files in `fastq_directory`. Reads are still filtered by `artic guppyplex` in both cases.
- Uncompressed FASTQ files (`.fastq` / `.fq`) are now accepted in the `fastq_1` and `fastq_2` columns for both Nanopore and Illumina samples.
- A warning is now raised when an Illumina sample has a `fastq_directory` value (which is ignored) or a Nanopore sample has a `fastq_2` value (which is ignored).
- A Nanopore sample which provides both `fastq_directory` and `fastq_1` now fails samplesheet validation with an informative error.

### `Fixed`

- Nanopore samples with an explicit FASTQ path and a `barcode` were silently dropped, explicit paths now always take precedence over `barcode` fuzzy matching.
- Fixed mismatched header / column order in the example `assets/samplesheet.csv`.

### `Dependencies`

### `Deprecated`
