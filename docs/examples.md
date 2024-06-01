# Examples

## Stripping very large JSON files

Some JSON files from Illumina TSO panels (for example, TSO500) are not QC
filtered and contain all detected genomic variants, irrespective of whether they
pass the quality criteria. Such files can get very large, too large to be
processed by any JSON parser. If your JSON file does not contain only
QC-filtered variants ("PASS"), it needs to be stripped (filtered) first before
using the icaparser module for further processing.

The code below can be run in Python in a terminal or in a Jupyter notebook.
Terminal is recommended.  

```python
import icaparser as icap
icap.strip_json_files(source_dir='../Data/Original', target_dir='../Data/Derived')
```

## Simple example

The code below is the *Hello World* example for reading and filtering ICA JSON
files with default filtering rules. For more sophisticated filtering options,
see the [API reference](reference.md).

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