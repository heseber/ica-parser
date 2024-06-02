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

## Installation of Illumina Connected Annotations

### Dowloading the software

Go to https://developer.illumina.com/illumina-connected-annotations, register and download the files. There is a ZIP file with binaries for Microsoft Dotnet, a Docker image, and some documentation. The Docker image is for arm64, i.e., it works on M1/M2 Macs, but not on Intel or AMD machines (Intel Macs, Linux).

### Installation alternatives

There are two options for installing _Illumina Connected Annotations_.

**[Option 1](Installing_ICA_nodocker.md):** Install the Microsoft Dotnet runtime and the _Illumina Connected Annotations_ binaries from the ZIP file

**[Option 2](Installing_ICA_docker.md):** Install the Docker image for your system and run _Illumina Connected Annotations_ from a Docker container.

## ICAparser

The Python package `icaparser` can be used to parse and filter the JSON files
that are the output of _Illumina Connected Annotations_.

The package provides functions for parsing JSON files created by Illumina's
Connected Annotations (ICA)pipeline.  ICA annotates mutations with a set of
tools and data sources. This package allows to:

* Strip JSON files from variants that do not pass quality criteria
* Load all mutations from a stripped Illumina Connected Annotations JSON file
* Filter mutations based on annotations and positions
* Create annotated tables of filtered mutations

See the help pages of the package for further details.

### Installation

#### Installation of the icaparser package

It is recommended to create a new virtual environment with Python >= 3.9 and to
install the icaparser package in that environment. Activate the environment and
run:

```sh
pip install "git+https://github.com/Bayer-Group/ica-parser.git"
```

If you want to install a particular development branch, use

```sh
pip install "git+https://github.com/Bayer-Group/ica-parser.git@BRANCHNAME"
```

If you use Jupyter notebooks, the virtual environment should be added as a new Jupyter kernel. See [Using Virtual Environments in Jupyter Notebook and Python - Parametric Thoughts](https://janakiev.com/blog/jupyter-virtual-envs/) how to do that.

#### Install ipywidgets

Required for progress bars in Jupyter. Please refer to the Jupyter or JupyterLab documentation how to install the widgets. For example:

```sh
conda install jupyter # if not installed yet
conda install jupyterlab_widgets
jupyter labextension install jupyter-matplotlib
jupyter lab build
exit
```

**→ Restart Jupyter**

### Usage

JSON files from Illumina's pipeline can be very large if they are not filtered for quality but contain any deviation from the reference genome, irrespective of the quality of the mutation call. Gzip compressed JSON files with sizes in the gigabyte range cannot be processed by JSON packages that read the entire file into memory. For such large volume JSON files it is necessary to first reduce the fle size by removing all variants that do not meet Illumina's quality criteria, i.e., by keeping only variants which have a `PASS` tag.

The stripped versions of the JSON files need to be created first before using any other function from this library.

### Example:

#### Step 1: Creating stripped JSON files

The code below can be run in Python in a terminal or in a Jupyter notebook. Terminal is recommended.  

```python
import icaparser as icap
icap.strip_json_files(source_dir='../Data/Original', target_dir='../Data/Derived')
```



#### Step 2: Processing stripped JSON files

The code below is the *Hello World* example for reading and filtering ICA JSON files with default filtering rules. For more sophisticated filtering options, see the help pages.

```python
import icaparser as icap
json_files = icap.get_dna_json_files('../Data/Derived')
first_file = json_files[0]
# Get the annotation data sources
icap.get_data_sources(first_file)
# Get pipeline run metadata
icap.get_pipeline_metadata(json_files)
# Get a mutation table
mut_table = icap.get_mutation_table_for_files(json_files)
```
