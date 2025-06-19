# artic-network/amplicon-nf pipeline parameters

Amplicon genome assembly which doesn\'t suck

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Type | Description | Default | Enum | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | `string` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>|  |  | True |  |
| `outdir` | `string` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. |  |  | True |  |
| `storedir` | `string` | Where to store files downloaded by the pipeline for future use. <details><summary>Help</summary><small>Any downloaded primer schemes / clair3 models will be stored in this directory for future use. If not set when running from the command line, the schemes will be stored in the current working directory/store_dir. If this is run via epi2me by default they will be stored in your `epi2melabs/data/` directory.</small></details>| ./storedir |  |  |  |
| `email` | `string` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>|  |  |  |  |
| `multiqc_title` | `string` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. |  |  |  |  |

## pipeline_settings

General settings which change the behaviour of the pipeline

| Parameter | Type | Description | Default | Enum | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| `normalise_depth` | `integer` | What depth of sequencing coverage to normalise the primertrimmed BAM files to. Set to 0 to disable normalisation. | 200 |  | True |  |
| `lower_ambiguity_frequency` | `number` | Minimum frequency of a minor allele variant to call an ambiguous variant. | 0.2 |  |  |  |
| `upper_ambiguity_frequency` | `number` | Maximum frequency of a minor allele variant to call an ambiguous variant. | 0.8 |  |  |  |
| `min_ambiguity_count` | `integer` | Minimum number of reads supporting a minor allele variant to call an ambiguous variant. | 1 |  |  |  |
| `min_coverage_depth` | `integer` | Minimum coverage depth for a position to be included in the consensus sequence. | 20 |  |  |  |
| `qc_pass_min_coverage` | `integer` | Coverage of reference above 'minimum_coverage_depth' to be considered a minimal QC pass | 50 |  |  |  |
| `qc_pass_high_coverage` | `integer` | Coverage of reference above `minimum_coverage_depth` to be considered a high QC pass. | 90 |  |  |  |
| `min_mapping_quality` | `integer` | Minimum mapping quality for a read to be included in the consensus sequence. | 20 |  |  |  |
| `primer_match_threshold` | `integer` | Allow fuzzy primer position matching within this threshold. | 35 |  |  |  |
| `min_ont_read_quality` | `integer` | Minimum mean quality score for a read to pass quality filtering. | 7 |  |  |  |
| `min_ont_read_length` | `integer` | Minimum length of a read to pass read length filtering. We recommend that you set this to your expected amplicon length - 100. |  |  |  |  |
| `max_ont_read_length` | `integer` | Maximum length of a read to pass read length filtering. We recommend that you set this to your expected amplicon length + 100. |  |  |  |  |
| `skip_ont_quality_check` | `boolean` | Skip filtering of ONT read lengths. <details><summary>Help</summary><small>Set this to true if you want to skip filtering of ONT read lengths. If you are providing barcode dirs from a 'fastq_pass' directory this should be set to true.</small></details>| False |  |  |  |
| `allow_mismatched_primers` | `boolean` | Do not remove reads which appear to have mismatched primers. | False |  |  |  |
| `manual_clair3_model` | `string` | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Clair3 model. <details><summary>Help</summary><small>Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/rerio/tree/master/clair3_models.</small></details>|  | ['r1041_e82_260bps_fast_g632', 'r1041_e82_400bps_fast_g632', 'r1041_e82_400bps_sup_g615', 'r1041_e82_260bps_hac_g632', 'r1041_e82_400bps_hac_g615', 'r1041_e82_400bps_sup_v400', 'r1041_e82_260bps_hac_v400', 'r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_sup_v410', 'r1041_e82_260bps_hac_v410', 'r1041_e82_400bps_hac_v400', 'r1041_e82_400bps_sup_v420', 'r1041_e82_260bps_sup_g632', 'r1041_e82_400bps_hac_v410', 'r1041_e82_400bps_sup_v430', 'r1041_e82_260bps_sup_v400', 'r1041_e82_400bps_hac_v420', 'r1041_e82_400bps_sup_v500', 'r1041_e82_260bps_sup_v410', 'r1041_e82_400bps_hac_v430', 'r104_e81_hac_g5015', 'r1041_e82_400bps_hac_v500', 'r104_e81_sup_g5015', 'r1041_e82_400bps_hac_v520', 'r1041_e82_400bps_sup_v520', 'r941_prom_sup_g5014', 'r941_prom_hac_g360+g422'] |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Type | Description | Default | Enum | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | `string` | Git commit id for Institutional configs. | master |  |  | True |
| `custom_config_base` | `string` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| https://raw.githubusercontent.com/nf-core/configs/master |  |  | True |
| `config_profile_name` | `string` | Institutional config name. |  |  |  | True |
| `config_profile_description` | `string` | Institutional config description. |  |  |  | True |
| `config_profile_contact` | `string` | Institutional config contact information. |  |  |  | True |
| `config_profile_url` | `string` | Institutional config URL link. |  |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Type | Description | Default | Enum | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| `version` | `boolean` | Display version and exit. |  |  |  | True |
| `publish_dir_mode` | `string` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| copy | ['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move'] |  | True |
| `email_on_fail` | `string` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>|  |  |  | True |
| `plaintext_email` | `boolean` | Send plain-text email instead of HTML. |  |  |  | True |
| `max_multiqc_email_size` | `string` | File size limit when attaching MultiQC reports to summary emails. | 25.MB |  |  | True |
| `monochrome_logs` | `boolean` | Do not use coloured log outputs. |  |  |  | True |
| `hook_url` | `string` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>|  |  |  | True |
| `multiqc_config` | `string` | Custom config file to supply to MultiQC. |  |  |  | True |
| `multiqc_logo` | `string` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file |  |  |  | True |
| `multiqc_methods_description` | `string` | Custom MultiQC yaml file containing HTML including a methods description. |  |  |  |  |
| `validate_params` | `boolean` | Boolean whether to validate parameters against the schema at runtime | True |  |  | True |
| `pipelines_testdata_base_path` | `string` | Base URL or local path to location of pipeline test dataset files | https://raw.githubusercontent.com/nf-core/test-datasets/ |  |  | True |
| `trace_report_suffix` | `string` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. |  |  |  | True |
