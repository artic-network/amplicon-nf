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

### Sample Reports

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>_amplicon-nf_sample-report.html`: a standalone sample report that can be viewed in your web browser, one report is generated for each sample provided in the samplesheet.

</details>

amplicon-nf generates custom per-sample reports with QC metrics and a read depth plot per genome segment (just one for non-segmented viruses) which is useful for investigating any issues encountered.

### Consensus FASTAs

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.consensus.fasta`: the produced consensus FASTA for each sample provided in the samplesheet.

</details>

The consensus FASTA contains the genome sequence of the sample based on the reads provided. This is based on the reference sequence provided, and areas not covered by sufficient reads to determine the contents of the genome are replaced with `N` indicating that any base may be present.

### BAM files

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.primertrimmed.sorted.bam`: The depth normalised, primertrimmed BAM file which was used to call variants against the reference.
  - `<SAMPLE>.sorted.bam`: Reads aligned to the reference with minimal filtering, normalisation, or primertrimming.

</details>

[BAM files](https://en.wikipedia.org/wiki/BAM_(file_format)) are an alignment format of reads aligned to one or multiple reference sequences. These can be useful to diagnose issues with your sequencing experiment.

### VCF files

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.vcf.gz`: The variant calls used to generate the consensus sequence in the FASTA.

</details>

[VCF files](https://en.wikipedia.org/wiki/Variant_Call_Format) are a text based format which describes how the consensus sequence is different from the reference sequence. The VCF file created will be slightly different for the Illumina / Nanopore workflows due to differences in how they are generated, most importantly the Illumina workflow calls mixed positions using [IUPAC ambiguity codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation) but since IUPAC codes are not valid VCF format a `ConsensusTag` `INFO` tag field indicates whether the position is `ambiguous` or `fixed`. For example:
```
MN908947.3      1875    .       C       T       .       .       DP=38;VAF=0.210526;ConsensusTag=ambiguous       .       .
```
Indicates that the position may be either `C` or `T` (IUPAC `Y`) not that a `C -> T` SNP was observed, the `VAF` field indicates the variant allele frequency observed.

### Amplicon depth TSV files

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/`
  - `<SAMPLE>.amplicon_depths.tsv`: The observed depths of each amplicon calculated by `align_trim`.

</details>

A [tab-separated-values](https://en.wikipedia.org/wiki/Tab-separated_values) file which summarises the mean depth of coverage for each amplicon in this sample. If no reads were assigned to an amplicon there will be no value in this TSV for that amplicon.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
