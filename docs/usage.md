# artic-network/amplicon-nf: Usage

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with at least 3 columns (more depending on how you are using the pipeline), and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

## Nanopore (ONT) only runs

There are two ways to run the pipeline for ONT data; 
1) with explicit FASTQ directories within the samplesheet in the `fastq_directory` column, this is recommended since all FASTQ directories are explicitly linked with their sample name within the samplesheet making mistakes where data is associated with incorrect metadata less likely.
2) with implicit (fuzzy) matching of FASTQ directories based on the provided `barcode` column, this option will match up subdirectories of the directory provided with the `--read_directory` parameter so you will not have to provide filepaths within the samplesheet.

> [!NOTE]
>Whichever option you use, it is important to remember that for ONT data the pipeline expects a directory of FASTQ files to be provided, not fastq files directly even if a directory only contains a single FASTQ file.

If you are only running ONT sequenced samples through the pipeline then you only need to fill in a subset of the samplesheet which will be different based on which of the above two options you are using, the following examples would be valid for a fastq_pass directory which looks like this:

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

### 1: Explicit FASTQ Directory Input

If you wish to provide explicit FASTQ directories then a valid samplesheet could look like this:
```csv title="samplesheet.nanopore_explicit.csv"
sample,platform,scheme_name,fastq_directory
sample1,nanopore,artic-inrb-mpox/2500/v1.0.0,/some/directory/fastq_pass/barcode01
sample2,nanopore,artic-inrb-mpox/2500/v1.0.0,/some/directory/fastq_pass/barcode02
sample3,nanopore,artic-inrb-mpox/2500/v1.0.0,/some/directory/fastq_pass/barcode03
```

> [!NOTE]
> An [example explicit Nanopore samplesheet](../assets/samplesheet.nanopore_explicit.csv) has been provided with the pipeline for reference.

### 2: Implicit (fuzzy) FASTQ Directory Input

If you wish to utilise fuzzy directory matching then a valid samplesheet could look like this (remember, the `fastq_pass` directory **MUST** be provided with `--read_directory` for this samplesheet to be valid):

```csv title="samplesheet.nanopore_fuzzy.csv"
sample,platform,scheme_name,barcode
sample1,nanopore,artic-inrb-mpox/2500/v1.0.0,barcode01
sample2,nanopore,artic-inrb-mpox/2500/v1.0.0,barcode02
sample3,nanopore,artic-inrb-mpox/2500/v1.0.0,barcode03
```

Please Note that the full barcode name must be provided and must match the directory exactly, as in, `01` or `1` would be invalid since they would not match the directory exactly.
> [!NOTE]
> An [example explicit Nanopore samplesheet](../assets/samplesheet.nanopore_implicit.csv) has been provided with the pipeline for reference.

## Illumina only runs

As with ONT data, there are two ways to run the pipeline for Illumina datasets; 
1) with explicit FASTQ directories within the samplesheet in the `fastq_1` and `fastq_2` columns, this is recommended since all FASTQ files are explicitly linked with their sample name within the samplesheet making mistakes where data is associated with incorrect metadata less likely.
2) with implicit (fuzzy) matching of FASTQ file pairs based on the provided `sample` column, this option will match up file pairs within the directory provided with the `--read_directory` parameter so you will not have to provide filepaths within the samplesheet.

If you are only running Illumina sequenced samples through the pipeline then you only need to fill in a subset of the samplesheet which will be different based on which of the above two options you are using, the following examples would be valid for a directory which looks like this:

```
run_fastq_directory
   ├── sample-1_S1_R1_001.fastq.gz
   |── sample-1_S1_R2_001.fastq.gz
   ├── sample-2_S2_R1_001.fastq.gz
   └── sample-2_S2_R2_001.fastq.gz
```
### 1: Explicit paired FASTQ input

For explicit FASTQ pairs, a valid samplesheet for the above directory would look like this:
```csv title="samplesheet.illumina_explicit.csv"
sample,platform,scheme_name,fastq_1,fastq_2
sample-1,illumina,artic-inrb-mpox/2500/v1.0.0,/some/directory/sample-1_S1_R1_001.fastq.gz,/some/directory/sample-1_S1_R2_001.fastq.gz
sample-2,illumina,artic-inrb-mpox/2500/v1.0.0,/some/directory/sample-2_S2_R1_001.fastq.gz,/some/directory/sample-2_S2_R2_001.fastq.gz
```

> [!NOTE]
> An [example Illumina samplesheet](../assets/samplesheet.illumina_explicit.csv) has been provided with the pipeline.

### 2: Fuzzy paired FASTQ Input

If you wish to utilise fuzzy directory matching then a valid samplesheet could look like this (remember, the `run_fastq_directory` **MUST** be provided with `--read_directory` for this samplesheet to be valid).

```csv title="samplesheet.illumina_implicit.csv"
sample,platform,scheme_name
sample-1,illumina,artic-inrb-mpox/2500/v1.0.0
sample-2,illumina,artic-inrb-mpox/2500/v1.0.0
```
> [!NOTE]
> An [example Illumina samplesheet](../assets/samplesheet.illumina_implicit.csv) has been provided with the pipeline.

### Primer Schemes

We recommend that you provide a scheme using a [primalscheme labs](https://labs.primalscheme.com/) identifier e.g. `artic-inrb-mpox/2500/v1.0.0` or `artic-sars-cov-2/400/v5.4.2` which is laid out with the following schema `<SCHEME_NAME>/<SCHEME_LENGTH>/<SCHEME_VERSION>`, the scheme itself will be sourced from the [primerschemes repository](https://github.com/quick-lab/primerschemes).

Alternatively, if you wish to use a scheme not available from the official repository you may provide a samplesheet containing the `custom_scheme_path` and `custom_scheme_name` parameters, `custom_scheme_path` should point to a directory containing two files `primer.bed` and `reference.fasta` which describe your custom scheme, `custom_scheme_name` is an optional field which allows you to provide a name for this custom scheme which will be used when generating a run report, if this is provided with a `scheme_name` the `custom_scheme_name` will be ignored. 

### Full samplesheet

If you wish to run a mix of platforms or even a mix of primer schemes in the same pipeline that is fully supported (even encouraged), if you provide a full samplesheet each sample is treated separately and a separate run QC report will be generated for each primer scheme.

A final samplesheet file consisting of both Nanopore and Illumina data and a mix of primer schemes may look something like the one below. 

```csv title="samplesheet.csv"
sample,platform,scheme_name,custom_scheme_path,custom_scheme_name,fastq_directory,fastq_1,fastq_2
nanopore_amplicon_data,nanopore,artic-inrb-mpox/2500/v1.0.0,,,/path/to/fastq/files/Barcode01/,,,
illumina_amplicon_data,illumina,,/path/to/custom_scheme/,some_scheme_name,,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
```

| Column               | Description                                                                                                                                                                                         |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`             | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`).              |
| `platform`           | The platform used to sequence the sample, this must either be `nanopore` or `illumina`                                                                                                              |
| `barcode`            | The barcode for this specific sample, if you provide this for nanopore data and the `--read_directory` parameter the pipeline will match your data with your samplesheet automatically              |
| `scheme_name`        | The [primalscheme labs](https://labs.primalscheme.com/) identifier for the scheme used to amplify your sequencing data, this must be in the format `<SCHEME_NAME>/<SCHEME_LENGTH>/<SCHEME_VERSION>` |
| `custom_scheme_path` | A path which points to a directory containing two files `primer.bed` and `reference.fasta` which describe your custom scheme.                                                                       |
| `custom_scheme_name` | The name of your custom scheme, this is used to refer to the scheme in the per-run QC reports.                                                                                                      |
| `fastq_directory`    | A directory containing your Nanopore read FASTQS for this sample.                                                                                                                                   |
| `fastq_1`            | Full path to FastQ file for Illumina short reads 1 for this sample.                                                                                                                                 |
| `fastq_2`            | Full path to FastQ file for Illumina short reads 2 for this sample.                                                                                                                                 |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run artic-network/amplicon-nf --input ./samplesheet.csv --outdir ./results --store_dir ./store_dir  -profile docker
```

This will launch the pipeline with the `docker` configuration profile which is recommended. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
<store_dir>          # Storage directory in specified location (defined with --store_dir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run artic-network/amplicon-nf -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
store_dir: './store_dir/'
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Running the pipeline in low resource settings

If you are running this pipeline locally (for example on a sequencing laptop) you may wish to put a limit on the amount of resources that the pipeline will attempt to use, to do this there are two profiles which limit the amount of resources the pipeline will try to use to fit on more modest hardware, if you wish to use one of these profiles you may do so like this:

```bash
nextflow run artic-network/amplicon-nf \
   -profile low_resource,<docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --store_dir <store_dir> 
```

The `-profile` parameter accepts multiple profiles separated by a comma so providing a parameter such as `-profile low_resource,docker` will use both profiles at the same time.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull artic-network/amplicon-nf
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [artic-network/amplicon-nf releases page](https://github.com/artic-network/amplicon-nf/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
