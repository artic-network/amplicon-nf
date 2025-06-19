# artic-network/amplicon-nf pipeline parameters

Amplicon genome assembly which doesn\'t suck

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `storedir` | Where to store files downloaded by the pipeline for future use. <details><summary>Help</summary><small>Any downloaded primer schemes / clair3 models will be stored in this directory for future use. If not set when running from the command line, the schemes will be stored in the current working directory/store_dir. If this is run via epi2me by default they will be stored in your `epi2melabs/data/` directory.</small></details>| `string` | ./storedir |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## pipeline_settings

General settings which change the behaviour of the pipeline

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `normalise_depth` | What depth of sequencing coverage to normalise the primertrimmed BAM files to. Set to 0 to disable normalisation. | `integer` | 200 | True |  |
| `lower_ambiguity_frequency` | Minimum frequency of a minor allele variant to call an ambiguous variant. | `number` | 0.2 |  |  |
| `upper_ambiguity_frequency` | Maximum frequency of a minor allele variant to call an ambiguous variant. | `number` | 0.8 |  |  |
| `min_ambiguity_count` | Minimum number of reads supporting a minor allele variant to call an ambiguous variant. | `integer` | 1 |  |  |
| `min_coverage_depth` | Minimum coverage depth for a position to be included in the consensus sequence. | `integer` | 20 |  |  |
| `qc_pass_min_coverage` | Coverage of reference above 'minimum_coverage_depth' to be considered a minimal QC pass | `integer` | 50 |  |  |
| `qc_pass_high_coverage` | Coverage of reference above `minimum_coverage_depth` to be considered a high QC pass. | `integer` | 90 |  |  |
| `min_mapping_quality` | Minimum mapping quality for a read to be included in the consensus sequence. | `integer` | 20 |  |  |
| `primer_match_threshold` | Allow fuzzy primer position matching within this threshold. | `integer` | 35 |  |  |
| `min_ont_read_quality` | Minimum mean quality score for a read to pass quality filtering. | `integer` | 7 |  |  |
| `min_ont_read_length` | Minimum length of a read to pass read length filtering. We recommend that you set this to your expected amplicon length - 100. | `integer` |  |  |  |
| `max_ont_read_length` | Maximum length of a read to pass read length filtering. We recommend that you set this to your expected amplicon length + 100. | `integer` |  |  |  |
| `skip_ont_quality_check` | Skip filtering of ONT read lengths. <details><summary>Help</summary><small>Set this to true if you want to skip filtering of ONT read lengths. If you are providing barcode dirs from a 'fastq_pass' directory this should be set to true.</small></details>| `boolean` | False |  |  |
| `allow_mismatched_primers` | Do not remove reads which appear to have mismatched primers. | `boolean` | False |  |  |
| `manual_clair3_model` | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Clair3 model. <details><summary>Help</summary><small>Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/rerio/tree/master/clair3_models.</small></details>| `string` |  |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  | True |
| `trace_report_suffix` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. | `string` |  |  | True |
