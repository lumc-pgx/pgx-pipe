# Pharmacogenomics analysis of long amplicon sequences
This Snakemake workflow assigns fully phased PGx haplotypes to targeted amplicon sequencing data generated using the PacBio platforms.


## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)
- [git](https://git-scm.com/)


## Installation
- Clone this repository
  - `git clone  https://github.com/lumc-pgx/pgx-pipe.git`

- Change to the pgx-pipe directory
  - `cd pgx-pipe`

- Initialize submodules
  - `git submodule init`

- Prepare submodules
  - `git submodule update`

- Create a conda environment for running the pipeline
  - `conda env create -n pxg-pipe -f environment.yaml`

- In order to use the pipeline on a computing cluster, update your .profile to use the drmaa library:
  - `echo "export DRMAA_LIBRARY_PATH=libdrmaa.so.1.0" >> ~/.profile`
  - `source ~/.profile`


## Configuration
The pipeline behavior is configured by editing [config.yaml](config.yaml).  


## Execution
- Activate the conda environment created during installation
  - `source activate pgx-pipe`

## Example:
  - `pipe-runner -s Snakefile -d output_dir -c config.yaml -cl cluster_settings.yaml`

Modify the config.yaml and cluster_settings.yaml according to the location of your data and run parameters.
The CYP2D6_locus_config.yaml is an example of locus definition file indicated in the config.yaml and should be changed according to the gene of interest.

- Note
  - The first run of the pipeline will create conda environments where required for the individual
    pipeline rules. This process can take some time.
  - Subsequent runs of the pipeline will re-use these environments, resulting in faster execution.

## Outputs
The output of the pipeline includes the following folders

```
- Annotation
- Haplotyping
- Phasing
- Structural variation
- Summary
- Variant calling

```


## Overview
Pgx-pipe combines the individual analysis modules into a meta-workflow which can be run with a single command.

          
