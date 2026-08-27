# Nextflow Pipeline for RNA-seq runs on CBSU BioHPC Servers 
[![](img/trex-extended-logo.png)](https://trex.biotech.cornell.edu/)

[![Build and Push Docker Image](https://github.com/bixBeta/nextflow/actions/workflows/docker.yml/badge.svg?branch=main)](https://github.com/bixBeta/nextflow/actions/workflows/docker.yml)
[![Docker Pulls](https://img.shields.io/docker/pulls/bixbeta/trex-rna)](https://hub.docker.com/r/bixbeta/trex-rna)
[![Docker Image Version](https://img.shields.io/docker/v/bixbeta/trex-rna/latest?label=docker%3Alatest)](https://hub.docker.com/r/bixbeta/trex-rna/tags)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![CZ ID](https://img.shields.io/badge/CZ%20ID-upload-6ab04c.svg)](https://czid.org/)

<hr>

This pipeline can be directly used on CBSU BioHPC Servers. To invoke the pipeline please make sure you have nextflow available in your path. 
You may use the [ following guide](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=638#c>) to help you configure nextflow for your account. 


To see if this pipeline works on your account simply run the following command on your terminal:

```
nextflow run bixBeta/nextflow -r main --help
```

If successful, you may see the following output on your terminal console:

```
R  N  A  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run bixBeta/nextflow -r main < args ... >

Args:
    * --listGenomes    : Get extended list of genomes available for this pipeline
    * --id             : TREx Project ID 
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -----------------------------------------------
        Sample Sheet Example: ( comma delimited file )    
        |-------|-----------------|-----------------|
        | label | fastq1          | fastq2          |
        |-------|-----------------|-----------------|
        | SS1   | SS1_R1.fastq.gz | SS1_R2.fastq.gz |
        |-------|-----------------|-----------------|
        | SS2   | SS2_R1.fastq.gz | SS2_R2.fastq.gz |
        |-------|-----------------|-----------------|
        .
        .
        . etc.
        -----------------------------------------------

    * --mode            : use 'PE'    for paired end data; default <PE>
                        : use 'PES'  for paired end data + split unmapped
                        : use 'PEB'  for paired end bacterial data
                        : use 'PEBS' for paired end bacterial data + split unmapped
                        : use 'SE'   for single end data
                        : use 'SES'  for single end data + split unmapped
                        : use 'SEB'  for single end bacterial data
                        : use 'SEBS' for single end bacterial data + split unmapped
                        : use 'SE3PL' for single end 3-prime poly-A library / Lexogen FWD (cutadapt + STAR)

    * --strand          : 0,1 or 2 for unstranded, first-strand and second-strand; default <2>
    * --fastqs          : Use this param if fastq files are in the fastqs folder in the project directory; 
                          If --fastqs is not specified, the fastqs must be supplied with absolute paths in the sample-sheet.csv 
    * --fastp           : Invokes fastp trimming module.
    * --genome          : Genome index. Use --listGenomes flag to see all available genomes. Also supports a path value for starIndex dir. 
    * --genome2         : Secondary Genome index. This will align the --genome subtracted reads to --genome2 index. (only use if --mode is PES, SES, PEBS or SEBS)
    * --screen          : Invokes the fastq_screen step. See the screen.conf file here <https://github.com/bixBeta/nextflow/blob/main/screen.conf> for more details. 
    * --gbcov           : Runs GeneBodyCoverage Program on sub-setted bams.
    * --chromosub       : Subset bams to specified chromosome name. < defaults to chromosome 10 >
    * --bed12           : Custom Path for BED12 file for geneBodyCoverage module
    * --splitname       : A string that will be used to denote --genome2 e.g. "GRC100011A", "Cat_custom" etc. 
    * --screenconf      : Supply custom screen config file, default ( <https://github.com/bixBeta/nextflow/blob/main/screen.conf> )
    * --mqcgenome       : A string denoting genome info. It will be used in the multiqc header and will also be used to organize results for a given run.
                          Highly recommended for genomes supplied as path to --genome param. e.g. --mqcgenome Cat_custom
    * --trinity         : Invokes Trinity de novo assembly on all fastp-trimmed reads (PE modes only); default <false>
                        : SS_lib_type is derived automatically from --strand (0=unstranded, 1=FR, 2=RF)
    * --salmon          : Quantify against Trinity assembly using Salmon (requires --trinity); default <false>
                        : Runs SuperTranscripts by default → gene-level quantification
    * --transcript_level : Use Trinity.fasta directly for Salmon (transcript-level instead of gene-level); requires --salmon
    * --czid            : Invokes CZ ID metagenomics upload subworkflow; default <false>
    * --czid_project    : CZ ID project name (must exist on czid.org)
    * --czid_host       : Host organism for CZ ID metadata; default <Homo sapiens>
                        : Accepts shorthand aliases: human, mouse, rat, chicken, dog, cat, cow, fly, mosquito, zebrafish, pig, rabbit, macaque
                        : Or pass any full CZ ID-recognised species name directly ( <https://czid.org/host_genomes> )
    * --czid_sample_type  : Sample type for CZ ID metadata; default <Tissue>
    * --czid_nucleotide   : Nucleotide type for CZ ID metadata; default <RNA> — use DNA for metagenomic DNA samples
```


<hr>
It is always a good idea to run the pull command before executing the pipeline. This will ensure that the user is using the latest branch version of the pipeline.
Use the following command to ensure the usage of the latest version of the pipeline:

```
nextflow pull bixBeta/nextflow -r main
```
<hr>

## Generate a Sample Sheet

Use `mk-sheet.py` to auto-generate a `sample-sheet.csv` from your Illumina delivery directory.
Labels are derived from `Sample_*` folder names — everything before the order number is kept (e.g. `Sample_8027D_BH10_10_104***_...` → `8027D_BH10_10`).

Download the script:

```bash
wget https://raw.githubusercontent.com/bixBeta/nextflow/main/mk-sheet.py
```

**Paired-end — absolute paths** (default):
```bash
python mk-sheet.py /local/Illumina/DRV/260***_RX_0***/Unaligned/Project_104***
```

**Paired-end — filenames only** (use alongside `--fastqs`, FASTQs must be in a `fastqs/` folder):
```bash
python mk-sheet.py /local/Illumina/DRV/.../Project_104*** --names
```

**Single-end:**
```bash
python mk-sheet.py /local/Illumina/DRV/.../Project_104*** --se
```

**Custom output filename:**
```bash
python mk-sheet.py /local/Illumina/DRV/.../Project_104*** -o my-sheet.csv
```

<hr>

## Parameters File (alternative approach)

Download the annotated `params.yaml` template to your working directory and edit it before running:

```bash
wget https://raw.githubusercontent.com/bixBeta/nextflow/main/params.yaml
```

Then run the pipeline using the params file instead of command-line flags:

```bash
nextflow run bixBeta/nextflow -r main -params-file params.yaml
```

> All parameters are documented with descriptions and defaults inside the file. Uncomment optional sections as needed.

<hr>
TREx Workflow:

![](img/trex-rna-tubemap.png)
