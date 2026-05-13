# Nextflow Pipeline for RNA-seq runs on GG02 
[![](img/trex-extended-logo.png)](https://trex.biotech.cornell.edu/)

[![Build and Push Docker Image](https://github.com/bixBeta/nextflow/actions/workflows/docker.yml/badge.svg?branch=czid)](https://github.com/bixBeta/nextflow/actions/workflows/docker.yml)
[![Docker Pulls](https://img.shields.io/docker/pulls/bixbeta/trex-rna)](https://hub.docker.com/r/bixbeta/trex-rna)
[![Docker Image Version](https://img.shields.io/docker/v/bixbeta/trex-rna/czid?label=docker%3Aczid)](https://hub.docker.com/r/bixbeta/trex-rna/tags)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![CZ ID](https://img.shields.io/badge/CZ%20ID-upload-6ab04c.svg)](https://czid.org/)

<hr>

This pipeline can be directly used on GG02. To envoke the pipeline please make sure you have nextflow available in your path. 
You may use the [ following guide](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=638#c>) to help you configure nextflow for your account. 


To see if this pipeline works on your account simply run the following command on your terminal:

```
nextflow run https://github.com/bixBeta/nextflow -r g2 --help
```

If successful, you may see the following output on your terminal console:

![](img/success.png)


<hr>
It is always a good idea to run the pull command before executing the pipeline. This will ensure that the user is using the latest branch version of the pipeline.
Use the following command to ensure the usage of the latest version of the pipeline:

```
nextflow pull https://github.com/bixBeta/nextflow -r g2 
```
<hr>
TREx Workflow:

![](img/trex-rna-tubemap.png)
