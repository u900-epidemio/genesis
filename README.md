# Variant calling pipeline 

**Variant calling pipeline for TUMOSPEC project**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.10-blue.svg)](https://multiqc.info/)


## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow manager to run tasks across multiple compute infrastructures in a very portable manner.
It supports [conda](https://docs.conda.io) package manager making installation easier and results highly reproducible.

## Pipeline summary

1. Run quality control of raw sequencing reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) 
2. Remove adapters of raw sequencing reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Align reads on reference genome ([`BWA-MEM`](https://github.com/lh3/bwa))
4. Filter reads on mapping quality ([`Samtools`](http://www.htslib.org/))
5. Intersection with the BED file ([`Bedtools`](https://bedtools.readthedocs.io/en/latest/))
6. Gene body coverage quality control ([`Samtools`](http://www.htslib.org/))
7. Variant calling ([`GATK`](https://gatk.broadinstitute.org/hc/en-us))
8. Variant filtering ([`GATK`](https://gatk.broadinstitute.org/hc/en-us))
9. VCF annotation with reference genome and databases: gnomAD(v2.1.1) and dbNSFP (v4.1a) ([`SnpEFF and SnpSIFT`](http://pcingola.github.io/SnpEff/))
10. VCF filtering with custom script ([`Python`](https://www.python.org/))
11. Present QC results in a final report ([`MultiQC`](https://multiqc.info/))



### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 19.10.0
Launching `main.nf` [stupefied_darwin] - revision: aa905ab621
=======================================================

Usage:
nextflow run main.nf --genome 'hg19' --samplePlan 'sample_plan_1.csv' -profile conda,cluster

Mandatory arguments:
--reads [file]                   Path to input data (must be surrounded with quotes)
--samplePlan [file]              Path to sample plan file if '--reads' is not specified
--genome [str]                   Name of the reference genome. See the `--genomeAnnotationPath` to defined the annotation path
-profile [str]                   Configuration profile to use (multiple profiles can be specified with comma separated values)
--snpeffDb [str]                 Directory to the snpEff databases

Skip options: All are false by default
--skipSoftVersion [bool]         Do not report software version
--skipMultiQC [bool]             Skip MultiQC

Other options:
--metadata [dir]                Add metadata file for multiQC report
--outDir [dir]                  The output directory where the results will be saved
-w/--work-dir [dir]             The temporary directory where intermediate data will be saved
-name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

=======================================================
Available profiles
-profile test                    Run the test dataset
-profile conda                   Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile multiconda              Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
-profile cluster                 Run the workflow on the cluster, instead of locally

```


### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follows:

#### Run the pipeline on a test dataset

See the file `conf/test.config` to set your test dataset. Download the snpEff database of the version declared as "nameSnpeff" in conf/genomes.conf in the "snpeffDb" directory. 

libcrypto.so might need to be changed to libcrypto.so.1.0.0 and exported in the LD_LIBRARY_PATH variable

```bash

nextflow run main.nf -profile test,conda --snpeffDb ${snpeffDb} -resume

```

#### Run the pipeline from a `sample plan`

```bash
nextflow run main.nf --samplePlan mySamplePlan.csv --genome 'hg38' --genomeAnnotationPath /my/annotation/path --outDir /my/output/dir --snpeffDb ${snpeffDb} -resume

```

### Defining the '-profile'

By default (whithout any profile), Nextflow excutes the pipeline locally, expecting that all tools are available from your `PATH` environment variable.

In addition, several Nextflow profiles are available that allow:
* the use of [conda](https://docs.conda.io) or containers instead of a local installation,
* the submission of the pipeline on a cluster instead of on a local architecture.

The description of each profile is available on the help message (see above).

Here are a few examples to set the profile options:

#### Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
```bash
-profile path --globalPath /my/path/to/bioinformatics/tools
```

#### Run the pipeline on the cluster, building a new conda environment
```bash
-profile cluster,conda --condaCacheDir /my/path/to/condaCacheDir

```

For details about the different profiles available, see [Profiles](docs/profiles.md).

### Sample plan

A sample plan is a csv file (comma separated) that lists all the samples with a biological IDs.
The sample plan is expected to contain the following fields (with no header):

```
SAMPLE_ID,SAMPLE_NAME,path/to/R1/fastq/file,path/to/R2/fastq/file (for paired-end only)
```


### Genome annotations

The pipeline does not provide any genomic annotations but expects them to be already available on your system. The path to the genomic annotations can be set with the `--genomeAnnotationPath` option as follows:

```bash
nextflow run main.nf --samplePlan mySamplePlan.csv --genome 'hg19' --genomeAnnotationPath /my/annotation/path --outDir /my/output/dir

```

For more details see  [Reference genomes](docs/referenceGenomes.md).

## Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/referenceGenomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie (S. Murat El Houdigui)

## Contacts

For any question, bug or suggestion, please use the issue system or contact the bioinformatics core facility.
