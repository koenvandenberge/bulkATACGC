# Normalization benchmark of bulk ATAC-seq data

[![DOI](https://zenodo.org/badge/269151763.svg)](https://zenodo.org/badge/latestdoi/269151763)


This repository contains code to reproduce analyses reported in our normalization benchmark paper.

To run the code, please first

- download the zipped file containing all public datasets from Zenodo at https://zenodo.org/record/6646500. The zipped file should be unzipped in the root directory of this repository, and the corresponding folder should be called 'data'. 
- download human GRCh38 genome at ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/. Download  mouse GRCm37.67 genome at ftp://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/ and GRCm38 genome at ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/. Set paths to these genomes in the R scripts accordingly.

The folders in this repository contain the following items.

- methods folder: contains scripts on custom normalization methods, or functions for the scone benchmark.
- objects folder: contains saved R objects. Either intermediate files that might be useful, or objects used for analysis.
- scripts folder: contains data analysis scripts.

Supplemental Data Figures are available from Zenodo at https://zenodo.org/record/6646500.

Additionally, simulated datasets for each of the real datasets may be downloaded from Zenodo at https://zenodo.org/record/6109369.
