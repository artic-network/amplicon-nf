<div align="center">
   <img src="/assets/amplicon-nf-badge.png" alt="artic-network/amplicon-nf" width="400">
</div>

[![GitHub Actions CI Status](https://github.com/artic-network/amplicon-nf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/artic-network/amplicon-nf/actions/workflows/nf-test.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17522200-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17522200)

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

## Acknowledgements

This pipeline has been created as part of the ARTIC network project funded by the Wellcome Trust (collaborator award – `313694/Z/24/Z` and discretionary award – `206298/Z/17/Z`) and is distributed as open source and open access. All non-code files are made available under a Creative Commons CC-BY licence unless otherwise specified. Please acknowledge or cite this repository or associated publications if used in derived work so we can provide our funders with evidence of impact in the field.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,barcode,platform,scheme_name,custom_scheme_path,custom_scheme_name,fastq_directory,fastq_1,fastq_2
nanopore_amplicon_data,,nanopore,artic-inrb-mpox/2500/v1.0.0,,,/path/to/fastq/files/Barcode01/,,,
illumina_amplicon_data,,illumina,,/path/to/custom_scheme/,some_scheme_name,,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
```

The `scheme_name` field refers to a scheme as a [primalscheme labs](https://labs.primalscheme.com/) identifier e.g. `artic-inrb-mpox/2500/v1.0.0` or `artic-sars-cov-2/400/v5.4.2`.

Each row represents a fastq file (single-end) or a pair of fastq files (paired end), the pipeline will run the Illumina and ONT workflows in parallel, it is important to note that the ONT and Illumina workflows have different input requirements. ONT requires only `fastq_directory` which is intended to be a directory as created by Dorado / minKNOW during basecalling, whereas Illumina requires a pair of read files.

> [!NOTE]
> There are more detailed pipeline / Nextflow usage instructions (including samplesheet construction and custom primer schemes), there are available in: [docs/usage.md](docs/usage.md).

Now, you can run the pipeline using:

```bash
nextflow run artic-network/amplicon-nf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --store_dir <STOREDIR> 
```

The pipeline is configured with a set of default parameters which should suit most use cases but a full list of available configurable parameters is available in [docs/parameters.md](https://github.com/artic-network/amplicon-nf/blob/main/docs/parameters.md).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Problems and Solutions

If you run into problems running this pipeline there is a list of known problems and their solutions [available in docs/problems.md](docs/problems.md). If the issue you encounter is not listed there please consider [making an issue!](https://github.com/artic-network/amplicon-nf/issues/new/choose)

## Credits

artic-network/amplicon-nf was originally written by Sam Wilkinson (@BioWilko).

I thank the following people for their assistance in the development of this pipeline:

* James A. Fellows Yates (@jfy133), Áine O'Toole (@aineniamh), Rachel Colquhoun (@rmcolq), Chris G. Kent (@ChrisgKent), Bede Constantinides (@bede) and Andrew Rambaut (@rambaut) for the extremely useful testing and feedback.
* Jared Simpson (@jts) for originally writing the [Illumina Freebayes consensus generation workflow](https://github.com/jts/ncov2019-artic-nf/blob/6ecf07bef462bfb896ae91629c49116761c03175/modules/illumina.nf#L174-L227) the Illumina workflow is based on.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use artic-network/amplicon-nf for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
