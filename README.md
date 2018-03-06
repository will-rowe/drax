<div align="center">
  <img src="assets/misc/drax-logo-with-text.png" alt="drax-logo" width="250">
  <h3><a style="color:#780200">D</a>etecting <a style="color:#780200">R</a>esistome <a style="color:#780200">A</a>ssociated ta<a style="color:#780200">X</a>a </h3>
  <hr>
  <a href="https://travis-ci.org/will-rowe/drax"><img src="https://travis-ci.org/will-rowe/drax.svg?branch=master" alt="travis"></a>
  <a href="https://www.nextflow.io"><img src="https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg" alt="Nextflow"></a>
  <a href="https://github.com/will-rowe/drax/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
</div>

***

`DRAX is currently under active development -- please check back soon.`

## Overview

`DRAX` is a [Nextflow](https://www.nextflow.io) pipeline to detect resistome associated taxa in metagenomes. It identifies Antimicrobial Resistance Genes (ARGs) using [GROOT](https://www.biorxiv.org/content/early/2018/02/28/270835), extracts the genomic environment of each ARG using [Metacherchant](https://academic.oup.com/bioinformatics/article-abstract/34/3/434/4575138?redirectedFrom=fulltext) and then performs taxanomic classification of assembled environments. Further pipeline steps and full documentation are being added.


## Quick Start

You can use [bioconda](https://bioconda.github.io/) to manage all of the `DRAX` pipeline dependencies:

```
# Get the software requirements list
wget https://raw.githubusercontent.com/will-rowe/drax/master/drax-conda-environment.yml

# Create a conda environment
yes | conda env create -n drax -f drax-conda-environment.yml && source activate drax

# Run DRAX
nextflow run will-rowe/drax --reads 'tests/*R{1,2}.fq.gz'
```

Alternatively, you can also run `DRAX` using containers:

```
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Run DRAX with Docker
./nextflow run will-rowe/drax --reads 'tests/*R{1,2}.fq.gz' -profile docker

# OR run DRAX with Singularity
./nextflow run will-rowe/drax --reads 'tests/*R{1,2}.fq.gz' -with-singularity 'docker://wpmr/drax'
```


## Documentation

The `DRAX` pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)
