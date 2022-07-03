# cas9-random-access
This repository contains the data processing and analysis pipeline for random access of DNA-archived data using Cas9 enrichment and nanopore sequencing.

## Included Packages
In the `software` directory, this repository includes an adapted version of the [C3POa software package](https://github.com/rvolden/C3POa) associated with the publication https://www.pnas.org/doi/10.1073/pnas.1806447115.

In `software/blat`, this repository includes a stand-alone version of the BLAT software from https://genome.ucsc.edu/cgi-bin/hgBlat.

## Dependencies
This code was developed and tested on Jupyter using Python 3.7.

The following packages should be installed separately:
- numpy (1.16.2)
- matplotlib (2.2.4)
- [poa (1.0.0 revision: 1.2.2.9)](https://github.com/tanghaibao/bio-pipeline)
- [gonk](https://github.com/rvolden/gonk)
- [minimap (2 2.7-r654)](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)

## Installation Instructions
Clone this repository and install the required dependencies.

In `software/config.txt`, specify the absolute path to the executables listed in the file.

In the second cell of the `20200715_triple_file_access.ipynb` notebook, replace `[PATH]` with the absolute path to the cloned repository.

## Run Example Analysis
