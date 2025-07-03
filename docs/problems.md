# artic-network/amplicon-nf: Problems & Solutions

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