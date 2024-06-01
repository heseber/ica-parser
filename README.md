# Parser for Illumina Connected Annotations JSON files

## Introduction

Illumina's _Illumina Connected Annotations_ (previous name: _Nirvana_) is a fast annotation tool for genomic alternations (VCF files) with similar functionality as the Variant Effect Predictor (VEP) or SnpEff.

Illumina decided not to continue releasing updates of their annotation pipeline as Open Source. The latest Open Source release is frozen to relatively old genome annotations (genes and transcripts) from 2017 and is not compatible with recent COSMIC versions. These annotations cannot be updated by users. Therefore, continued use of the Open Source version is discouraged.

Illumina released binary-only versions of their internal pipeline. This pipeline uses up-to-date genome annotations and is now also compatible with recent COSMIC releases. Compared to the last public release of Nirvana (v3.18.1) this release includes new features such as: 

- Updated Transcripts from RefSeq and Ensembl
- New annotation sources:
  - Primate AI-3D
  - Cancer Hotspots
- Updated annotations (Oct 2023)
  - ClinVar
  - OMIM
  - ClinGen
  - dbSNP 156
- MANE select tags and Consequence impacts for Transcript annotations
- More accurate HGVS notations

PrimateAI-3D, spliceAI, COSMIC, and OMIM annotations are available from Illumina by request only because they require a license for commercial use. 

## Dowloading the software

Go to https://developer.illumina.com/illumina-connected-annotations, register and download the files. There is a ZIP file with binaries for Microsoft Dotnet, a Docker image, and some documentation. The Docker image is for arm64, i.e., it works on M1/M2 Macs, but not on Intel or AMD machines (Intel Macs, Linux).

## Installation alternatives

There are two options for installing _Illumina Connected Annotations_.

**[Option 1](Installing_ICA_nodocker.md):** Install the Microsoft Dotnet runtime and the _Illumina Connected Annotations_ binaries from the ZIP file

**[Option 2](Installing_ICA_docker.md):** Install the Docker image for your system and run _Illumina Connected Annotations_ from a Docker container.

## ICAparser

The Python package [icaparser](package/README.md) can be used to parse and filter the JSON files that are the output of _Illumina Connected Annotations_.
