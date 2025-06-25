# artic-network/amplicon-nf

[![GitHub Actions CI Status](https://github.com/artic-network/amplicon-nf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/artic-network/amplicon-nf/actions/workflows/nf-test.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/artic-network/amplicon-nf)

## Introduction

**artic-network/amplicon-nf** is a bioinformatics pipeline that takes sequencing reads generated from ARTIC-style viral amplicon sequencing schemes, assembles them into consensus sequences, and runs some basic quality control on the outputs.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->



<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,platform,scheme_name,custom_scheme_path,custom_scheme_name,fastq_directory,fastq_1,fastq_2
nanopore_amplicon_data,nanopore,artic-inrb-mpox/2500/v1.0.0,,,/path/to/fastq/files/Barcode01/,,,
illumina_amplicon_data,illumina,,/path/to/custom_scheme/,some_scheme_name,,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end), the pipeline will run the Illumina and ONT workflows in parallel, it is important to note that the ONT and Illumina workflows have different input requirements. ONT requires only `fastq_directory` which is intended to be a directory as created by Dorado / minKNOW during basecalling. Below there is an example layout of a `fastq_pass` directory, each row of the samplesheet in this case would point to a single `barcode` directory.

```
fastq_pass
   ├── barcode01
   |   ├── reads0.fastq.gz
   │   └── reads1.fastq.gz
   ├── barcode02
   │   ├── reads0.fastq.gz
   │   ├── reads1.fastq.gz
   │   └── reads2.fastq.gz
   └── barcode03
       └── reads0.fastq.gz
```

For example, if your `fastq_pass` directory looked the same as the above example your samplesheet could look like this:

```csv
sample,platform,scheme_name,fastq_directory
barcode01,nanopore,artic-inrb-mpox/2500/v1.0.0,/some/directory/fastq_pass/barcode01
barcode02,nanopore,artic-inrb-mpox/2500/v1.0.0,/some/directory/fastq_pass/barcode02
barcode03,nanopore,artic-inrb-mpox/2500/v1.0.0,/some/directory/fastq_pass/barcode03
```

Please note that for a single platform run (as in, your data is all ONT rather than a mix of ONT and Illumina) you do not need to include irrelevant columns, the same is true for the custom scheme fields (`custom_scheme_path` and `custom_scheme_name`).

For Illumina data the paths to the forward and reverse read files should be provided directly, for example if your directory looks like the below:
```
run_fastq_directory
   ├── sample-1_S1_L001_R1_001.fastq.gz
   |── sample-1_S1_L001_R2_001.fastq.gz
   ├── sample-2_S2_L001_R1_001.fastq.gz
   └── sample-2_S2_L001_R2_001.fastq.gz
```

Then your samplesheet could look like this:

```csv
sample,platform,scheme_name,fastq_1,fastq_2
sample-1,illumina,artic-inrb-mpox/2500/v1.0.0,/some/directory/sample-1_S1_L001_R1_001.fastq.gz,/some/directory/sample-1_S1_L001_R2_001.fastq.gz
sample-2,illumina,artic-inrb-mpox/2500/v1.0.0,/some/directory/sample-2_S2_L001_R1_001.fastq.gz,/some/directory/sample-2_S2_L001_R2_001.fastq.gz
```

As the ONT only samplesheet above does not need the `fastq_1` or `fastq_2` path columns the Illumina only samplesheet does not need the `fastq_directory` column.

We recommend that you provide a scheme using a [primalscheme labs](https://labs.primalscheme.com/) identifier e.g. `artic-inrb-mpox/2500/v1.0.0` or `artic-sars-cov-2/400/v5.4.2` which is laid out with the following schema `<SCHEME_NAME>/<SCHEME_LENGTH>/<SCHEME_VERSION>`, the scheme itself will be sourced from the [primerschemes repository](https://github.com/quick-lab/primerschemes).

Alternatively, if you wish to use a scheme not available from the official repository you may provide a samplesheet containing the `custom_scheme_path` and `custom_scheme_name` parameters, `custom_scheme_path` should point to a directory containing two files `primer.bed` and `reference.fasta` which describe your custom scheme, `custom_scheme_name` is an optional field which allows you to provide a name for this custom scheme which will be used when generating a run report, if this is provided with a `scheme_name` the `custom_scheme_name` will be ignored. 

Now, you can run the pipeline using:

```bash
nextflow run artic-network/amplicon-nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --storedir <STOREDIR> 
```

The pipeline is configured with a set of default parameters which should suit most use cases but a full list of available configurable parameters is available in [docs/parameters.md](https://github.com/artic-network/amplicon-nf/blob/main/docs/parameters.md).

If you are running this pipeline locally (for example on a sequencing laptop) you may wish to put a limit on the amount of resources that the pipeline will attempt to use, to do this there are two profiles which limit the amount of resources the pipeline will try to use to fit on more modest hardware, if you wish to use one of these profiles you may do so like this:

```bash
nextflow run artic-network/amplicon-nf \
   -profile low_resource,<docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --storedir <STOREDIR> 
```

The `-profile` parameter accepts multiple profiles separated by a comma so providing a parameter such as `-profile low_resource,docker` will use both profiles at the same time.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Problems and Solutions

### No Basecall Model in Read Headers

If you provide ONT data which does not have the `basecall_model_version_id` field in the read header (ONT basecallers include this automatically) you will get an error message that looks like this:

```sh
  Provided fastq does not contain basecall_model_version_id in the read header so clair3 model cannot be chosen automatically, please provide an appropriate model with the --model parameter
```

If you do see this, you will need to provide the Clair3 model name manually with the `--manual_clair3_model` parameter, a full list of the available models is available in the [parameters.md document](docs/parameters.md), please note that this will apply to all samples in the run so a single samplesheet should only include data generated using a single basecalling model if using the `--manual_clair3_model` parameter.

### Insufficient Memory / CPUs available

All processes in this pipeline have a set amount of `cpus` and `memory` it requests, we have tried to make these requests as reasonable as possible while remaining realistic but in some circumstances you might see an error message like this:
```sh
Caused by:
  Process requirement exceeds available memory -- req: 8 GB; avail: 7.8 GB
```
If you do, we have included two profiles which will cap how much memory or CPUs a process can request `-profile low_resource`  and `-profile lower_resource`, these profiles will cap memory requests to `15Gb` and `7Gb` respectively, they will also reduce the number of processes which can run in parallel to attempt to reduce the overall load on your system.

These profiles may still lead to crashes if you have a large run to process, this is unfortunately unavoidable, if you run into such errors then we recommend that you try splitting your samplesheet into smaller sub runs.

## Credits

artic-network/amplicon-nf was originally written by Sam Wilkinson.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use artic-network/amplicon-nf for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
