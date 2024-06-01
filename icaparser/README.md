# Parser for Illumina Connected Annotations JSON files

This Python package provides functions for parsing JSON files created by Illumina's Connected Annotations (ICA)pipeline.  ICA annotates mutations with a set of tools and data sources. This package allows to:

* Strip JSON files from variants that do not pass quality criteria
* Load all mutations from a stripped Illumina Connected Annotations JSON file
* Filter mutations based on annotations and positions
* Create annotated tables of filtered mutations

See the help pages of the package for further details.

## Installation

### Installation of the icaparser package

It is recommended to create a new virtual environment with Python >= 3.9 and to install the icaparser package in that environment. Activate the environment and run:

```sh
pip install "git+https://github.com/Bayer-Group/ica-parser.git#subdirectory=icaparser"
```

If you want to install a particular development branch, use

```sh
pip install "git+https://github.com/Bayer-Group/ica-parser.git@BRANCHNAME#subdirectory=icaparser"
```

If you use Jupyter notebooks, the virtual environment should be added as a new Jupyter kernel. See [Using Virtual Environments in Jupyter Notebook and Python - Parametric Thoughts](https://janakiev.com/blog/jupyter-virtual-envs/) how to do that.

## Install ipywidgets

Required for progress bars in Jupyter. Please refer to the Jupyter or JupyterLab documentation how to install the widgets. For example:

```sh
conda install jupyter # if not installed yet
conda install jupyterlab_widgets
jupyter labextension install jupyter-matplotlib
jupyter lab build
exit
```

**→ Restart Jupyter**

## Usage

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
