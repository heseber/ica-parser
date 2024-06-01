import os

import pandas as pd
from icaparser.icaparser import (
    strip_json_file,  # noqa: F401
    strip_json_files,  # noqa: F401
    get_dna_json_files,  # noqa: F401
    get_header,  # noqa: F401
    get_sample,  # noqa: F401
    get_header_scalars,  # noqa: F401
    get_data_sources,  # noqa: F401
    get_pipeline_metadata,  # noqa: F401
    get_positions,  # noqa: F401
    get_multi_sample_positions,  # noqa: F401
    get_genes,  # noqa: F401
    get_position_by_coordinates,  # noqa: F401
    get_max_af,  # noqa: F401
    get_gnomad_max_af,  # noqa: F401
    get_gnomad_exome_max_af,  # noqa: F401
    get_onekg_max_af,  # noqa: F401
    get_cosmic_max_sample_count,  # noqa: F401
    get_clinvar_max_significance,  # noqa: F401
    cleanup_cosmic,  # noqa: F401
    common_variant_filter,  # noqa: F401
    filter_positions_by_variants,  # noqa: F401
    filter_variants_by_transcripts,  # noqa: F401
    filter_positions_by_transcripts,  # noqa: F401
    get_clinvar,  # noqa: F401
    get_consequences,  # noqa: F401
    get_vep_rank_for_consequence,  # noqa: F401
    get_vep_priority_for_consequence,  # noqa: F401
    get_vep_consequence_for_rank,  # noqa: F401
    get_strongest_vep_consequence_rank,  # noqa: F401
    get_strongest_vep_consequence_name,  # noqa: F401
    get_strongest_vep_consequence_priority,  # noqa: F401
    get_biotype_priority,  # noqa: F401
    get_gene_type,  # noqa: F401
    get_mutation_table_for_position,  # noqa: F401
    get_mutation_table_for_positions,  # noqa: F401
    get_mutation_table_for_files,  # noqa: F401
    explode_consequence,  # noqa: F401
    get_default_gene_type_map,  # noqa: F401
    get_default_mutation_classification_rules,  # noqa: F401
    apply_mutation_classification_rules,  # noqa: F401
    get_default_mutation_aggregation_rules,  # noqa: F401
    get_aggregated_mutation_table,  # noqa: F401
    split_multi_sample_json_file,  # noqa: F401
)

__version__ = "0.2.4"
__author__ = "Henrik Seidel"
