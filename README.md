# cas9-random-access
This repository contains the data processing and analysis pipeline for random access of DNA-archived data using Cas9 enrichment and nanopore sequencing.

## Included Packages
Our data processing pipeline uses an adapted version of the [C3POa software package](https://github.com/rvolden/C3POa) associated with the publication https://www.pnas.org/doi/10.1073/pnas.1806447115. This is found in the `software` directory.

In `software/blat`, a stand-alone version of the BLAT software from https://genome.ucsc.edu/cgi-bin/hgBlat is included. BLAT is called by C3POa to perform alignment. 

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

In the second cell of the `20200715_triple_file_access.ipynb` notebook, replace `[PATH]` with the absolute path to this cloned repository.

## Example Analysis
The directory `20200715_basecalled_reads` contains FASTQ files of nanopore reads that passed basecalling. These reads came from an experiment where 3 files (file #2, #13, and #24) were accessed from a pool of 25 files.

An example analysis of these reads is presented in the notebook `20200715_triple_file_access.ipynb`. In the Preprocessing and Processing sections of this notebook, the reads are demultiplexed into their respective files and a consensus sequence is called for each concatemeric read using C3POa. These consensus sequences are stored in FASTA format in `20200715_results/consensus`. The Analysis section of this notebook uses these consensus reads to generate relevant figures.
