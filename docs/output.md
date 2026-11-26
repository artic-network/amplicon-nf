# artic-network/amplicon-nf: Output

<!-- ## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. -->

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

<!-- ## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

-

- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution -->

### Run Reports

<details markdown="1">
<summary>Output files</summary>

- `<SCHEME_NAME>_amplicon-nf_run-report.html`: a standalone run report that can be viewed in your web browser, one report is generated per primer scheme used in the pipeline.

</details>

amplicon-nf generates custom run reports for each primer scheme used to generate the data used in the run. This report summarises some basic QC information about each sample in a table as well as some plots which are helpful to determine the cause of any issues encountered.

The QC table shows two coverage metrics.
- **Raw Coverage** is the percentage of the whole reference with depth ≥ `--min_coverage_depth` before normalisation, and is what the QC result is based on.
- **Amplicon Coverage** is the same measure restricted to the regions amplicons can cover: the union of amplicon inserts, with primers excluded. It uses post-normalisation depth when normalising.

Any run-level warnings are shown in a red box at the top of the report. For example, a warning is raised when pre- and post-normalisation amplicon coverage disagree, which usually means `--normalise_depth` is lower than `--min_coverage_depth`.

The **Primer Trimming QC** tab summarises the per-sample `align_trim` QC reports (see [Amplicon depth and primer trimming QC files](#amplicon-depth-and-primer-trimming-qc-files)) across all samples of the scheme. Values are shown as reported, with no pass/fail thresholds. The tab has five sub-tabs:
- **Read Fates**: kept and removed fragments by reason.
- **Amplicons**: per-amplicon fragment fates.
- **Primer Pairs**: primer pair types.
- **Split Reads**: split read join types.
- **Overlaps**: overlap ambiguity and base-frequency comparison.

### Run QC Summaries

<details markdown="1">
<summary>Output files</summary>

- `<SCHEME_NAME>_amplicon-nf_qc_results.tsv`: a TSV file containing the information in the run report QC table plus the primer scheme and the full pre/post-normalisation breakdown. That breakdown covers whole-genome coverage, amplicon coverage (`pre_norm_amplicon_coverage` / `post_norm_amplicon_coverage`) and mean depth.

</details>

A TSV ([tab-separated-values](https://en.wikipedia.org/wiki/Tab-separated_values)) file containing basic summary QC info about the samples included in the run, a summary TSV is generated per scheme used in the run.

### Sample Reports

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>_amplicon-nf_sample-report.html`: a standalone sample report that can be viewed in your web browser, one report is generated for each sample provided in the samplesheet.
  - `<SAMPLE>_amplicon-nf_qc-report.tsv`: a TSV file containing the per-contig QC summary for the sample, with the same columns as the run QC summary TSV.

</details>

amplicon-nf generates custom per-sample reports with QC metrics and a read depth plot per genome segment (just one for non-segmented viruses) which is useful for investigating any issues encountered. A companion TSV is also written for each sample containing the QC summary table (`sample`, `contig`, `contig_alias`, `primer_scheme`, `pre_norm_coverage`, `pre_norm_mean_depth`, `post_norm_coverage`, `post_norm_mean_depth`, `total_reads`, `total_amp_dropouts`, `amplicon_dropouts`, `qc_result`), with one row per genome segment. `amplicon_dropouts` is a comma-separated list of the amplicons with a mean depth below `--min_coverage_depth`.

### Consensus FASTAs

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.consensus.fasta`: the produced consensus FASTA for each sample provided in the samplesheet.

</details>

The consensus FASTA contains the genome sequence of the sample based on the reads provided. This is based on the reference sequence provided, and areas not covered by sufficient reads to determine the contents of the genome are replaced with `N` indicating that any base may be present.

### Combined Consensus FASTAs

<details markdown="1">
<summary>Output files</summary>

- `<CHROM>.<SCHEME_NAME>.combined_consensus.fasta`: A consensus FASTA for each segment of the sequenced virus (just one for non-segmented viruses) for each scheme included in the run.

</details>

The combined consensus FASTA contains the genome sequence of all samples for a specific genome segment based on the reads provided. This is based on the reference sequence provided, and areas not covered by sufficient reads to determine the contents of the genome are replaced with `N` indicating that any base may be present.

### Aligned FASTAs

<details markdown="1">
<summary>Output files</summary>

- `<CHROM-NAME>.<SCHEME>.aligned-consensus.fasta`: An alignment of all consensus sequences aligned to the reference sequence for each reference segment for each scheme.

</details>

The alignment FASTA(s) contain(s) an alignment of all consensus FASTAs to the scheme reference FASTA with one FASTA file per segment and per scheme.

### BAM files

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.primertrimmed.sorted.bam`: The depth normalised, primertrimmed BAM file which was used to call variants against the reference.
  - `<SAMPLE>.sorted.bam`: Reads aligned to the reference with minimal filtering, normalisation, or primertrimming.

</details>

[BAM files](<https://en.wikipedia.org/wiki/BAM_(file_format)>) are an alignment format of reads aligned to one or multiple reference sequences. These can be useful to diagnose issues with your sequencing experiment.

### VCF files

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.vcf.gz`: The variant calls used to generate the consensus sequence in the FASTA.
  - `<SAMPLE>.fail.vcf` (Nanopore) / `<SAMPLE>.fail.vcf.gz` (Illumina): Variants which were called but failed filtering, so their position is masked as `N` in the consensus sequence.
    - Nanopore: the `artic minion` failed variants, written as Clair3 reported them (the `FILTER` column is not updated). A variant fails if its quality is below `--min_variant_quality`, its allele frequency is between `--min_mask_allele_frequency` and `--min_allele_frequency`, it is a low quality frameshift, or its depth is below `--min_coverage_depth`.
    - Illumina: the mixed positions (`ConsensusTag=ambiguous`, see below), which are also included in `<SAMPLE>.vcf.gz`. With `--allow_iupac_codes` these positions are represented by IUPAC ambiguity codes rather than masked.

</details>

[VCF files](https://en.wikipedia.org/wiki/Variant_Call_Format) are a text based format which describes how the consensus sequence is different from the reference sequence. The VCF file created will be slightly different for the Illumina / Nanopore workflows due to differences in how they are generated, most importantly the Illumina workflow records mixed positions (variants with an allele frequency between `--min_mask_allele_frequency` and `--min_allele_frequency`) in the VCF. In both workflows these positions are masked as `N` in the consensus sequence by default, but the Illumina workflow can instead represent them using [IUPAC ambiguity codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation) if `--allow_iupac_codes` is set. Since IUPAC codes are not valid VCF format a `ConsensusTag` `INFO` tag field indicates whether the position is `ambiguous` or `fixed`. For example:

```
MN908947.3      1875    .       C       T       .       .       DP=38;VAF=0.210526;ConsensusTag=ambiguous       .       .
```

Indicates that the position may be either `C` or `T` (IUPAC `Y`) not that a `C -> T` SNP was observed, the `VAF` field indicates the variant allele frequency observed.

### Read QC and trimming (Illumina only)

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.fastp.json`: [fastp](https://github.com/OpenGene/fastp) report in JSON format, summarising read quality and adapter trimming statistics before/after trimming.
  - `<SAMPLE>.fastp.html`: Standalone HTML report of the same fastp statistics.
  - `<SAMPLE>.fastp.log`: fastp execution log.
  - `fastqc/`
    - `<SAMPLE>_raw_fastqc.html` / `.zip`: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) report for the raw, untrimmed reads.
    - `<SAMPLE>_trim_fastqc.html` / `.zip`: FastQC report for the reads after fastp trimming.

</details>

Illumina reads are quality- and adapter-trimmed with [fastp](https://github.com/OpenGene/fastp), which auto-detects adapter content rather than requiring a known adapter sequence. FastQC is run both before and after trimming so that the effect of trimming on read quality and adapter content can be assessed. All of these reports are summarised together in the [MultiQC](#multiqc) report.

### Reference IUPAC masking (Illumina only)

Some primer scheme reference sequences contain [IUPAC ambiguity codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation) (e.g. `R`, `Y`, `W`, `S`, `K`, `M`, `B`, `D`, `H`, `V`) to indicate positions where the reference itself is not a single definitive base. The variant caller used for Illumina samples ([Freebayes](https://github.com/freebayes/freebayes)) does not support IUPAC codes in the reference and will fail if they are present.

When the pipeline detects IUPAC ambiguity codes in a reference sequence it automatically replaces them with `N` before alignment and variant calling, and prints the following warning to the pipeline log:

```
WARN: IUPAC ambiguity codes were detected in one or more reference sequences and have been replaced with N...
```

This masking affects alignment sensitivity at those positions — reads covering a masked base are slightly less likely to map with high confidence because `N` matches any base equally — but this is preferable to incorrect variant calls or a pipeline failure. If your reference contains a large number of IUPAC positions you may wish to resolve them to definitive bases manually before running the pipeline.

### Amplicon depth and primer trimming QC files

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.amplicon_depths.tsv`: The observed depths of each amplicon calculated by `align_trim`.
  - `<SAMPLE>.summary.tsv`: The total number of read fragments, the number kept, and the number removed for each removal reason (e.g. `low_mapq`, `no_primer_match`, `ambiguous_amplicon`, `mismatched_primers`, `by_normalise`).
  - `<SAMPLE>.amplicons.tsv`: The number of fragments assigned, kept and removed (per removal reason) for each amplicon in the scheme.
  - `<SAMPLE>.primer_pairs.tsv`: Counts of each combination of left / right primers matched at the fragment ends, classified by type (e.g. `same_amplicon`, `spans_amplicons_same_pool`, `inside_overlap`).
  - `<SAMPLE>.split_reads.tsv`: Split (supplementary) alignments and where the parts join, classified by type (`amplicon_chimera`, `joined_inside_amplicon`, `low_mapq_part`).
  - `<SAMPLE>.overlaps.tsv`: For each overlap between neighbouring amplicons, the fragments removed as ambiguous and the base frequency concordance between the two amplicons covering it.

</details>

The amplicon depths file is a [tab-separated-values](https://en.wikipedia.org/wiki/Tab-separated_values) file which summarises the mean depth of coverage for each amplicon in this sample. If no reads were assigned to an amplicon there will be no value in this TSV for that amplicon.

The remaining TSV files are QC reports written by `align_trim` during primer trimming and normalisation, for both Nanopore and Illumina samples. They describe what happened to every read, which is useful when investigating amplicon dropouts, low coverage or chimeric reads. They are summarised across samples in the run report's Primer Trimming QC tab.

### Nextclade per-sample outputs

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.nextclade_report.tsv`: Nextclade QC metrics and clade/lineage assignment for the sample in tabular format.
  - `<SAMPLE>.nextclade_report.json`: Full Nextclade output in JSON format, including detailed mutation and QC information.

</details>

Per-sample Nextclade results are written into each sample's output directory alongside the other per-sample files (BAM, VCF, consensus FASTA, etc.).

A table summarising Nextclade QC results across all samples is also included in the run report (`<SCHEME_NAME>_amplicon-nf_run-report.html`).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pathoplexus submission folders

<details markdown="1">
<summary>Output files</summary>

- `pathoplexus/<SCHEME_NAME>/`: one folder per primer scheme (`/` in the scheme name is replaced with `_`). Not written with `--skip_pathoplexus_template`.
  - `<SCHEME_NAME>_pathoplexus_metadata.xlsx`: the metadata template. It has three sheets:
    - `Data`: one row per sample. Highlighted columns are for you to fill in, and auto-filled values are in grey italics.
    - `_lists`: hidden, holds the dropdown values.
    - `Guidance`: instructions, field definitions and provenance.
  - `<SCHEME_NAME>_pathoplexus_metadata.tsv`: the same rows and columns as the `Data` sheet. Use this for API submissions.
  - `<SCHEME_NAME>_pathoplexus_sequences.fasta`: the consensus sequences, renamed to match the metadata.
    - A sample with one sequence is named `<SAMPLE>`.
    - A segmented sample has one record per segment, named `<SAMPLE>_<REFERENCE>`, and these are listed in the `fastaIds` column.
    - Entirely-`N` sequences are left out.
  - `raw_reads/`: the dehosted reads of every sample in the template, named as in `files.rawReads`, with FASTQ headers stripped to the read ID. Not written with `--skip_dehosting`.

</details>

These folders hold what you need to submit a scheme's samples to [Pathoplexus](https://pathoplexus.org). See [Pathoplexus submissions](pathoplexus.md) for the fields you need to fill in, how every pre-filled field is derived and the scope of this feature: it is **not** designed for primer panels or multi-pathogen schemes, where the reference contains several sequences that are not segments of one genome.

- **Reference-selecting schemes:** artic maps each sample to the closest variant reference for some schemes (e.g. `artic-pan-ebola`). The primers are the same, so the folder and `ampliconPcrPrimerScheme` use the scheme name as given in the samplesheet.
- **Disk space:** `raw_reads/` holds header-stripped copies of the dehosted reads, so they take up disk space again.
- **Cloud executors:** the reads are copied into and out of a task, roughly three times the read volume in transfer.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.
  - Clade/Lineages are output as a tsv: `lineages.tsv`

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
