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

- For parallel execution on the cluster, writing output to the default *output* directory
  - `pipe-runner`

- To specify that the pipeline should write output to a location other than the default:
  - `pipe-runner --directory path/to/output/directory`

- Note
  - The first run of the pipeline will create conda environments where required for the individual
    pipeline rules. This process can take some time.
  - Subsequent runs of the pipeline will re-use these environments, resulting in faster execution.


## Overview
Pgx-pipe combines the individual analysis modules into a meta-workflow which can be run with a single command.

          
