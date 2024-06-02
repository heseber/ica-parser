# Welcome to the icaparser Python package

This Python package provides functions for parsing JSON files created by
Illumina's Connected Annotations (ICA)pipeline.  ICA annotates mutations with ‚a
set of tools and data sources. This package allows to:

* Strip JSON files from variants that do not pass quality criteria
* Load all mutations from a stripped Illumina Connected Annotations JSON file
* Filter mutations based on annotations and positions
* Aggregate mutations to gene level
* Create annotated tables of filtered mutations

See the [examples](examples.md) and the [API documentation](reference.md) for
further details.


