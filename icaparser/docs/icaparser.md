Module icaparser
====================
Parser for JSON files from Illumina's ICA annotation pipeline.

Functions
---------

    
`add_gene_types(positions)`
:   Adds the gene type to each transcript.
    
    Transcripts will be annotated with the gene type (*oncogene*, *tsg*,
    *mixed*) by adding a new attribute `geneType`. Only transcripts with
    one of these three gene types get this additional annotation. Other
    transcripts will not get the `geneType` attribute.
    
    Args:
        positions: list of filtered or unfiltered positions from JSON files
    
    Returns:
        A list of positions with additional annotation of transcripts
    
    Examples:
        >>> import icaparser as icap
        >>> positions = icap.add_gene_types(positions)

    
`apply_mutation_classification_rules(positions, rule_set={'gof': {'mutated': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}, 'uncertain': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}}, 'lof': {'mutated': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}, 'uncertain': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}}}, gene_type_map={'oncogene': {'gof'}, 'tsg': {'lof'}, 'mixed': {'gof', 'lof'}, 'Act': {'gof'}, 'LoF': {'lof'}, 'ambiguous': {'gof', 'lof'}, '': {'gof', 'lof'}}, hide_progress=False)`
:   Applies mutation classification rules to all positions.
    
    Each mutation is categorized for each isoform that overlaps with the genomic
    position of the mutation. Each transcript which passed the "mutated" or
    "uncertain" rule of the classification rules gets a new attribute
    `mutation_status` with the value "mutated" or "uncertain".
    
    In addition to classifying mutations based on their impact on transcript
    isoforms, this function also assembles the mutation status on sample and
    gene level without further aggregation. The impact depends on the type of
    gene ("gof" or "lof"), so the impacts are assembled separately for each gene
    type:
    
    `sample_id` → `gene` → `gene_type` → `variant_id` → `mutationStatus`
    
    The impact of a particular mutational variant can be different for different
    overlapping transcript variants of a gene, and the transcript variants can
    also belong to different genes. The strongest impact on any overlapping
    transcript of a gene is defined as the impact of that mutational
    variant on the gene. The analyst must decide which isoforms are used to
    classify genes. For example, only canonical transcripts may be considered.
    Alternatively, all transcripts or a subset of transcripts may be used.
    Therefore, it is necessary to first apply transcript-level filters to all
    genomic positions  before determining the mutation status of genes.
    
    Args:
        positions: list of positions
        rule_set: rules for classifying "gof" and "lof" genes. See also
            `get_default_mutation_classification_rules()` for an example
            if the default rule set needs to be modified.
        gene_type_map: dictionary for mapping gene types to canonical gene
            types. See also `get_default_gene_type_map()` for an example
            if the default mapping table needs to be modified.
    
    Returns:
        A list of positions and a dictionary with assembled and aggregated
        mutations
    
    Examples:
        >>> import icaparser as icap
        >>> positions, sample_muts = icap.apply_mutation_classification_rules(positions)

    
`cleanup_cosmic(positions)`
:   Remove Cosmic entries with alleles not matching the variant alleles.
    
    ICA attaches Cosmic entries to variants based on position only, which
    leads to wrong assignments of Cosmic entries to variants. This function
    removes all Cosmic entries from a variant for which reference and altered
    alleles do not match those of the variant.
    
    Filtering is done in place.
    
    Args:
        positions: list of positions to clean up
    
    Returns:
        A list of positions with cleaned up Cosmic entries

    
`common_variant_filter(variant, max_af=0.01)`
:   Get a variant filter based on GnomAD, GnomAd Exome, and 1000 Genomes.
    
    Returns True if none of the maximum allele frequencies from GnomAD,
    GnomAD Exome and 1000 Genomes is larger than `max_af`.
    
    Args:
        variant: the variant to investigate
        max_af: the mixum allele frequences threshold
    
    Returns:
        A bool

    
`explode_consequence(mutation_table, inplace=False)`
:   Explode the VEP consequence column of a mutation table.
    
    Exploding the VEP consequence column with the standard Pandas
    `explode()` function would return consquences as strings, not as
    ordered categories. This function will instead return a consequence
    columns which is an ordered category. The categories are ordered by
    their impact.
    
    Args:
        mutation_table: the mutation table to explode
        inplace: if True, then modify the mutation_table in place instead
            of returning a new object
    
    Returns:
        A pandas.DataFrame
    
    Examples:
        >>> import icaparser as icap
        >>> icap.explode_consequence(mutation_table, inplace=True)
    
        >>> mutation_table_exploded = icap.explode_consequence(mutation_table)

    
`filter_positions_by_transcripts(positions, filter_func)`
:   Filter positions based on a filter function for transcripts.
    
    Apply a filter function to all transcripts of each position.
    Transcripts not passing the filter are removed from a position.
    Positions without any transcripts passing the filter are removed
    from the returned list.
    
    Args:
        positions: list of positions to filter
        filter_func: function taking a transcript and returning a bool.
                     True means to keep the transcript.
    
    Returns:
        A list
    
    Example:
        >>> is_canonical_transcript = lambda x: x.get('isCanonical', False)
        >>> canonical_positions = icap.filter_positions_by_transcripts(
                non_common_positions,
                is_canonical_transcript
            )

    
`filter_positions_by_variants(positions, filter_func)`
:   Filter positions based on a filter function for variants.
    
    Apply a filter function to all variants of each position.
    Variants not passing the filter are removed from a position.
    Positions without any variants passing the filter are removed
    from the returned list.
    
    Args:
        positions: list of positions to filter
        filter_func: function taking a variant and returning a bool.
                     True means to keep the variant.
    
    Returns:
        A list
    
    Example:
        >>> import icaparser as icap
        >>> ax_af = 0.01
        >>> is_not_common_variant = lambda x: icap.common_variant_filter(x, max_af)
        >>> non_common_positions = icap.filter_positions_by_variants(
                positions,
                is_not_common_variant
            )

    
`filter_variants_by_transcripts(variants, filter_func)`
:   Filter variants based on a filter function for transcripts.
    
    Apply a filter function to all transcripts of each variant.
    Transcripts not passing the filter are removed from a variant.
    Variants without any transcripts passing the filter are removed
    from the returned list.
    
    Args:
        variants: list of variants to filter
        filter_func: function taking a transcript and returning a bool.
                     True means to keep the transcript.
    
    Returns:
        A list

    
`get_aggregated_mutation_table(positions, mutation_classification_rules={'gof': {'mutated': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}, 'uncertain': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}}, 'lof': {'mutated': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}, 'uncertain': {'position_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'variant_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>, 'transcript_filter': <function get_default_mutation_classification_rules.<locals>.<lambda>>}}}, mutation_aggregation_rules={'gof': <function get_default_mutation_aggregation_rules.<locals>._aggregation_rules_gof>, 'lof': <function get_default_mutation_aggregation_rules.<locals>._aggregation_rules_lof>}, gene_type_map={'oncogene': {'gof'}, 'tsg': {'lof'}, 'mixed': {'gof', 'lof'}, 'Act': {'gof'}, 'LoF': {'lof'}, 'ambiguous': {'gof', 'lof'}, '': {'gof', 'lof'}}, hide_progress=False)`
:   Returns a sample-gene-mutationStatus table.
    
    This function applies mutation classification rules to all mutational
    variants and aggregates the mutations according to the aggregation rules.
    This results in a table with one row for each sample-gene pair. The table
    contains several columns with impacts for "lof" and "gof" on allele level
    and gene level and with one additional column with the maximum impact for
    both allele and gene level.
    
    Args:
        positions: list of positions
        mutation_classification_rules: rules for classifying single mutations.
            See `get_default_mutation_classification_rules()` for details.
        mutation_aggregation_rules: rules for aggregation mutations.
            See `get_default_mutation_aggregation_reles()` for details.
        gene_type_map: dictonary for mapping gene types to canonical gene types.
            See `get_default_gene_type_map()` for details.
    
    Returns:
        A pandas.DataFrame with the mutation table.

    
`get_biotype_priority(biotype)`
:   Get the numeric priority of a biotype.
    
    The numeric priority of a biotype that is returned by this function
    is the same as defined by https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl.
    Biotypes are 'protein_coding', 'LRG_gene', ,'miRNA', ...
    
    Args:
        biotype: a string with a single biotype
    
    Returns:
        An integer, smaller values mean higher priority

    
`get_clinvar(variant)`
:   Get a table of all ClinVar annotations for a variant.
    
    Args:
        variant: the variant to investigate
    
    Returns:
        A pandas.DataFrame

    
`get_clinvar_max_significance(variant, ordered_significances=['pathogenic', 'likely pathogenic', 'drug response', 'risk factor', 'protective', 'affects', 'association', 'uncertain significance', 'conflicting data from submitters', 'other', 'not provided', 'likely benign', 'benign'])`
:   Get the maximum signifinance for all ClinVar annotations of a variant.
    
    Args:
        variant: the variant to investigate
        ordered_significances: ranked order of ClinVar significances
    
    Returns:
        A string

    
`get_consequences(transcript)`
:   Get a list of consequences for a transcript
    
    A list of consequences for a transcript is returned. If any of the annotated
    consequences is a combination if single consequences, separated by
    ampersands ('&') or commas, such a consequence is split into single consequences.
    
    Args:
        transcript: the transcript to get the consequences for
    
    Returns:
        A list of strings

    
`get_cosmic_cancer_gene_census_tier(variant)`
:   Get the highest Cosmic Cancer Gene Census tier of a variant.
    
    A variant can have no, one or multiple associated Cosmic identifiers.
    This function returns the highest 'tier' for all of the Cosmic Cancer
    Gene Census entries for a variant. There are two tiers, '1', and '2'.
    Tier '1' is assigned to genes with documented activity relevant to cancer,
    and tier '2' is assigned to genes with strong indications of a role in
    cancer but with less extensive available evidence.
    
    It is important to understand that this is the tier for a gene, not for a
    particular variant, so it must not be mistaken with the Cosmic Cancer
    Mutation Census tiers.
    
    See https://cancer.sanger.ac.uk/census for a detailed description of tiers.
    
    Args:
        variant:    the variant to investigate
    
    Returns:
        An integer

    
`get_cosmic_max_sample_count(variant, only_fully_annotated=False, only_allele_specific=True)`
:   Get the maximum sample count for all Cosmic annotations of a variant.
    
    A variant can have no, one or multiple associated Cosmic identifiers.
    This function returns the maximum sample count of all Cosmic identifers.
    For each Cosmic identifiers, sample numbers are summed up across all indications.
    Returns 0 if no Cosmic identifier exists for this variant.
    
    The 'only_fully_annotated' argument can be used to exclude Cosmic entries
    that have no 'cancerTypesAndCounts' annotation. If, for example, the
    Cosmic VCF file used by ICA contained all Cosmic variants, but the
    Cosmic TSV file used by ICA contained only Cancer Census variants, then
    a variant that is in the VCF but not the TSV will have a simply annotation
    without 'cancerTypesAndCounts', 'cancerSitesAndCounts', 'tiersAndCounts'.
    So when setting 'only_fully_annotated' to True, only samples having a mutation
    from the Cancer Census will be counted.
    
    The 'only_allele_specific' argument is used to exclude Cosmic entries that
    annotate the same chromosomal location but an allele that is different from
    the allele of the annotated variant. ICA annotates a variant with all
    Cosmic entries for that chromosomal location, irrespective of alleles.
    When counting Cosmic samples, this leads to an overestimation of Cosmic
    sample counts for a particular variant. Therefore, 'only_allele_specific'
    is True by default to count only samples from Cosmic entries with matching
    alleles. Occasionally, it may be desired, though, to count all samples with
    mutations at a given position, irrespective of allele. For example, several
    different alleles at a functional site of a gene can lead to function-disrupting    
    mutations, so we want to get the maximum sample count for any allele at
    that position. One might also think of adding the samnple counts for all
    Cosmic entries annotating a variant, but this does not work currently due
    to redundancy of Cosmic entries. Older Cosmic versions often included the
    same sample in different Cosmic entries. And newer Cosmic versions often
    have multiple entries for an allele, one for each transcript variant, with
    the same underlying samples.
    
    Args:
        variant:              the variant to investigate
        only_fully_annotated: consider only consmic entries which are fully
                              annotated, i.e. which have 'cancerTypesAndCounts'
        only_allele_specific: consider only cosmic entries with alleles
                              matching the allele of the annotated variant
    
    Returns:
        An integer

    
`get_data_sources(file)`
:   Extract a table with annotation data sources from the JSON header.
    
    Args:
        file: name of the ICA JSON file
    
    Returns:
        A pandas.DataFrame

    
`get_default_gene_type_map()`
:   Returns the default gene type map.
    
    The canonical gene types are `gof`, `lof`, and the union of both.
    Genes that need to be activated to drive a tumor are of type `gof`.
    Genes that need to be deactivated to drive a tumor are of type `lof`.
    Genes that need to be activated or deactivated depending on the context
    are of the union of both types.
    Genes for which it is unknown if they need to be activated or deactivated
    are also annotated with both types.
    Genes can be originally annotated with other type names than the canonical
    ones. The gene type map is used to map these other gene type names to the
    canonical gene types.
    
    The default map is:
    
    - `oncogene` → `{"gof"}`
    - `tsg` → `{"lof"}`
    - `Act` → `{"gof"}`
    - `LoF` → `{"lof"}`
    - `mixed` → `{"gof", "lof"}`
    - `ambiguous` → `{"gof", "lof"}`
    
    Returns:
        A dictionary with mappings from gene types to canonical gene types
    
    Examples:
        >>> import icaparser as icap
        >>> icap.get_default_gene_type_map()

    
`get_default_mutation_aggregation_rules()`
:   Returns the default mutation aggregation rules.
    
    Two types of the mutation status of a gene are introduced - allele level and
    gene level:
    
    - For "gof" genes (like oncogenes) it is sufficient if one of the alleles of
        one of the relevant isoforms has an activating mutation.
    - For "lof" genes (like tumor suppressor genes) all alleles of all relevant
        isoforms need to be functionally disrupted, either by mutations or by
        other means.
    - For Mixed and Other genes, the impact of a mutation is defined as the
        highest impact according to "gof" rules and "lof" rules.
    
    For "gof" genes, the allele and gene level classifications are identical
    unless there is additional information about activating modifications other
    than mutations. For "lof" genes, allele and gene level classifications
    may be different. For example, a truncating mutation of a tumor suppressor
    gene is functionally disruptive for the affected allele. However, there may
    be other functionally active alleles of the same gene, i.e., the gene itself
    may still be active. All alleles of such a gene must be dysfunctional,
    either by additional mutations or by other processes, such as copy number
    deletions or hypermethylation. Therefore, a single variant that is
    disruptive at the allele level is not necessarily also disruptive at the
    gene level. For "lof" genes, we usually do not have enough information
    for a reliable estimation of functional effects. Instead, some heuristic
    ules must be applied, and the analyst must decide whether to work with
    allele-level or gene-level classifications. We designate a "lof" gene as
    functionally disrupted (strong impact) if it has at least two mutations with
    either strong impact or uncertain impact. If a "lof" gene has only one of
    these mutations, it is designated as uncertain impact at the gene level,
    even if one of these mutations has a strong impact at the allele level. By
    categorizing effects at both allele and gene levels, we retain the
    flexibility to decide in downstream analyses how to merge some of these
    categories for subsequent statistical calculations.
    
    Additional filters based on publications, white lists and black lists for
    mutations can be applied. These are not part of the first version and will
    be added later. White and black lists can be based, for example, on ClinVar
    and on the publication by Hess et al. (2019).
    
    Returns:
        A dictionary
    
    Examples:
        >>> import icaparser as icap
        >>> icap.get_default_mutation_aggregation_rules()

    
`get_default_mutation_classification_rules(cosmic_threshold=10)`
:   Returns the default rules for classifying mutations.
    
    Defines the default rules for classifying mutations. The returned dictonary
    has keys "gof" and "lof", and the respective values are the rule sets for
    these gene types. A rule set is again a dictionary with the keys "mutated"
    and "uncertain". The values for "mutated" or "uncertain" are dictionaries
    with three filter functions, a "position_filter", a "variant_filter", and
    a "transcript filter". For example, a transcript will be called "mutated"
    if all three filters for "mutated" return True, and it will be called
    "uncertain", if all three filter functions for "uncertain" return True.
    
    Args:
        cosmic_threshold: for "gof" genes, this is the "hotspot threshold" for
            Cosmic, i.e., the minimum number of samples in Cosmic having that
            mutation to consider a mutation a hot spot and, therefore, call the
            mutation "mutated". If the number of Cosmic samples is smaller, the
            mutation is called "uncertain".
    
    Returns:
        A dictionary with mutation classification rules
    
    Examples:
        >>> import icaparser as icap
        >>> icap.get_default_mutation_classification_rules()
        >>> icap.get_default_mutation_classification_rules(cosmic_threshold=20)

    
`get_dna_json_files(base_dir, pattern='*MergedVariants_Annotated_filtered.json.gz')`
:   Find DNA annotation JSON files in or below `base_dir`.
    
    Searches for ICA DNA annotation JSON files in and below `base_dir`.
    All file names matching `pattern` are returned.
    
    Args:
        base_dir: base directory of directory subtree where to search
                  for DNA annotation JSON files
        pattern:  files names matching this pattern are returned
    
    Returns:
        A list of file names

    
`get_gene_type(gene_symbol)`
:   Get the gene type (oncogene, tsg, mixed) for a gene.
    
    Args:
        gene_symbol: the gene symbol of the gene
    
    Returns:
        a string

    
`get_genes(file)`
:   Extract gene annotation from a ICA JSON file.
    
    The `genes` section of ICA JSON files is optional. If this section
    is not included in the file, an empty list is returned.
    
    Args:
        file: name of the ICA JSON file
    
    Returns:
        A list

    
`get_gnomad_exome_max_af(variant, cohorts=['afr', 'amr', 'eas', 'nfe', 'sas'])`
:   Get the maximum allele frequency for gnomAD Exome.
    
    Get the maximum allele frequences across all major cohorts annotated
    by gnomAD, Exome  excluding bottleneck populations (Ashkenazy Jews and
    Finish) and _other_.
    
    Args:
        variant: the variant to investigate
        cohorts: subpopulations to include
    
    Returns:
        A float

    
`get_gnomad_max_af(variant, cohorts=['afr', 'amr', 'eas', 'nfe', 'sas'])`
:   Get the maximum allele frequency for gnomAD.
    
    Get the maximum allele frequences across all major cohorts annotated
    by gnomAD, excluding bottleneck populations (Ashkenazy Jews and Finish)
    and _other_.
    
    Args:
        variant: the variant to investigate
        cohorts: subpopulations to include
    
    Returns:
        A float

    
`get_header(file)`
:   Extract the header element from a ICA JSON file.
    
    Args:
        file: name of the ICA JSON file
    
    Returns:
        A dictionary with the header from the JSON file

    
`get_header_scalars(file)`
:   Extract a table with all scalar attributes from the JSON header.
    
    Args:
        file: name of the ICA JSON file
    
    Returns:
        A pandas.DataFrame

    
`get_max_af(variant, source, cohorts=None)`
:   Get the maximum allele frequency for a particular annotation source.
    
    Get the maximum allele frequency across all cohorts annotated by the
    annotation source.
    
    Args:
        variant: the variant to investigate
        source: the annotation source to use
        cohorts: subpopulations to include; include all if cohorts == None
    
    Returns:
        A float
    
    Example:
        >>> import icaparser as icap
        >>> icap.get_max_af(variant, 'gnomad')

    
`get_multi_sample_positions(files, *args, **kwargs)`
:   Extract all positions for a set of ICA JSON files.
    
    The sample id is stored as an additional new attribute of the
    `samples` element of a position. The `samples` element is a list,
    although ICA usually only creates single sample JSON files.
    
    Args:
        files: names of the ICA JSON files
        args: extra arguments forwarded to get_positions()
        kwargs: extra named arguments forwarded to get_positions()
    
    Returns:
        A list
    
    Examples:
        >>> import icaparser as icap
        >>> positions = icap.get_multi_sample_positions(json_files)
        >>> print(positions[0]['samples'][0]['sampleId'])

    
`get_mutation_table_for_files(json_files, max_af=0.01, min_vep_consequence_priority=6, min_cosmic_sample_count=0, only_canonical=False, extra_variant_filters=[], extra_transcript_filters=[])`
:   Get an annotated table of all filtered transcripts from a list of ICA JSON files.
    
    Load all positions from a list of ICA JSON files and filter them.
    Positions having any remaining variants and transcripts passing the filter
    are returned as an annotated table.
    
    Args:
        files:  list of ICA JSON files
        max_af: maximum allele frequency for gnomAD, gnomAD Exome and 1000 Genomes.
                Only variants with maximum allele frequencies below this threshold
                will be returned.
        min_vep_consequence_priority: only transcripts with a minimum VEP
                consequence priority not larger than this threshold will be retained
        min_cosmic_sample_count: only variants with a maximum cosmic sample count
                not lower than this threshold will be retained
        only_canonical: if true, only canonical transcripts will be retained
        extra_variant_filters: any additional filters to apply to variants.
                Filters shall return True to keep a variant.
        extra_transcript_filters: any additional filters to apply to transcripts.
                Filters shall return True to keep a transcript.
    
    Returns:
        A pandas.DataFrame
    
    Examples:
        >>> import icaparser as icap
        >>> extra_transcript_filters = [
                lambda x: x.get('source', '') == 'Ensembl',
                lambda x: x.get('hgnc', '') == 'KRAS'
            ]
        >>> mut_table = icap.get_mutation_table_for_files(
                json_files,
                extra_transcript_filters=extra_transcript_filters
            )

    
`get_mutation_table_for_position(position)`
:   Get an annotated table of all transcripts.
    
    Returns an annotated table of all transcripts that are affected
    by a mutation at a position.
    
    Args:
        position: the position to investigate
    
    Returns:
        A pandas.DataFrame

    
`get_mutation_table_for_positions(positions, hide_progress=False)`
:   Get an annotated table of all transcripts for all positions.
    
    Returns an annotated table of all transcripts that are affected
    by a mutation at any of the position.
    
    Args:
        positions: the positions to investigate
    
    Returns:
        A pandas.DataFrame

    
`get_onekg_max_af(variant)`
:   Get the maximum allele frequency for the 1000 Genomes Project.
    
    Get the maximum allele frequences across all cohorts annotated
    by the 1000 Genomes Project.
    
    Args:
        variant: the variant to investigate
    
    Returns:
        A float

    
`get_pipeline_metadata(files)`
:   Extract a table with metadata annotation pipeline run from the JSON header.
    
    Args:
        file: name of the ICA JSON file
    
    Returns:
        A pandas.DataFrame

    
`get_position_by_coordinates(positions, chromosome, position)`
:   Extract a particular position from a position list.
    
    Args:
        positions: a list of positions
        chromosome: name of the chromosome
        position: numeric posion on the chromosome
    
    Returns:
        A dict
    
    Examples:
        >>> import icaparser as icap
        >>> icap.get_position_by_coordinates(positions, 'chr1', 204399064)

    
`get_positions(file, variant_filters=[], transcript_filters=[])`
:   Extract all positions from a ICA JSON file.
    
    The sample id is stored as an additional new attribute of the
    `samples` element of a position. The `samples` element is a list,
    although ICA usually only creates single sample JSON files.
    
    Args:
        file: name of the ICA JSON file
        variant_filters: any filters to apply to variants.
                Filters shall return True to keep a variant.
        transcript_filters: any filters to apply to transcripts.
                Filters shall return True to keep a transcript.
    
    Returns:
        A list
    
    Examples:
        >>> transcript_filters = [
                lambda x: x.get('source', '') == 'Ensembl',
                lambda x: x.get('hgnc', '') == 'KRAS'
            ]
        >>> positions = icap.get_sample_positions(
                json_file,
                transcript_filters = transcript_filters
            )
        >>> print(positions[0]['samples'][0]['sampleId'])

    
`get_sample(file, suffix='(-D[^.]*)?\\.bam')`
:   Extract the sample name from a ICA JSON file.
    
    Args:
        file: name of the ICA JSON file
        suffix: regular expression to remove from the sample name in the JSON
            file. Defaults to '(-D[^.]*)?\.bam'.
    
    Returns:
        A string with the name of the sample annotated in the JSON file

    
`get_strongest_vep_consequence_name(transcript)`
:   Get the name of the strongest VEP consequence for a transcript.
    
    Args:
        transcript: the transcript to investigate
    
    Returns:
        A string

    
`get_strongest_vep_consequence_priority(transcript)`
:   Get the strongest priority of VEP consequence for a transcript.
    
    Get the strongest numeric priority of all VEP consequences
    for a transcript. Smaller numeric priorities mean stronger impact.
    
    Args:
        transcript: the transcript to investigate
    
    Returns:
        An integer

    
`get_strongest_vep_consequence_rank(transcript)`
:   Get the strongest rank of VEP consequences for a transcript.
    
    Get the strongest numeric rank of all VEP consequences
    for a transcript. Smaller ranks mean stronger impact.
    
    The priority of consequences is taken into account first. So if
    two consequences have different priorities, the consequence with
    the higher priority (lower priority number) will be used, and the
    rank for this consequence will be returned. If there are multiple
    consequences with the same priority, the lowest (strongest) rank
    will be returned.
    
    For clarification: ranks are unique, i.e. all VEP consequences
    ordered as listed on the VEP documentation page get the row number
    of this table assigned as rank.
    
    However, several consequences can have the same priority (e.g., stop
    gained and frameshift have the same priority). Priorities are copied
    from vcf2maf.pl.
    
    Args:
        transcript: the transcript to investigate
    
    Returns:
        An integer

    
`get_vep_consequence_for_rank(rank)`
:   Get the VEP consequence term of a numeric rank.
    
    Args:
        rank: the numeric rank of the consequence term
    
    Returns:
        A string

    
`get_vep_priority_for_consequence(consequence)`
:   Get the numeric priority of a VEP consequence term.
    
      The numeric priority of a consequence that is returned by this function
      is the same as defined by https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl.
    
    Args:
          consequence: the consequence term of the mutation
    
      Returns:
          An integer

    
`get_vep_rank_for_consequence(consequence)`
:   Get the numeric rank of a VEP consequence term.
    
    The numeric rank of a consequence that is the position of the consequence
    in this list of consequences at https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    
    Args:
        consequence: the consequence term of the mutation
    
    Returns:
        An integer

    
`split_multi_sample_json_file(json_file, output_dir)`
:   Splits a multi-sample JSON file into sample specifc JSON files.
    
    This function reads a multi-sample JSON file that was generated by
    annotating a multi-sample VCF file with ICA and splits it into
    sample-specific JSON files.
    
    Annotating very many single-sample VCF files with ICA is very time
    consuming, because ICA reads all annotation sources for each VCF file
    and this is dominating the runtime of ICA. It is therefore helpful to
    first merge many single-sample VCF files into one or a small number of
    multi-sample VCF files (for example, with `bcftools merge`), to annotate
    the multi-sample VCF file with ICA, and then to split the multi-sample
    JSON output of ICA into single-sample JSON files. These single-sample
    JSON files are required for the rest of this package.
    
    Args:
        json_file: the multi-sample json input file
        output_dir: the directory where to write the single sample JSON files.
            The directory will be created if it does not exist.
    
    Returns:
        Nothing

    
`strip_json_file(ifname, ofname)`
:   Reduce the JSON file size by keeping only 'PASS' variants.
    
    JSON files from Illumina's ICA pipeline can be very large because they
    contain any deviation from the reference genome, irrespective of the
    quality of the mutation call. Gzip compressed JSON files with sizes in the
    gigabyte range cannot be processed by JSON packages that read the entire
    file into memory. It is necessary to first reduce the size of JSON files by
    removing all variants that do not meet Illumina's quality criteria.
    
    This function reads a single JSON file and creates a single JSON outpout
    file by removing all variants that do not pass Illumina's quality criteria.
    
    Args:
        ifname: name of the  input file
        ofname: name of the output file
    
    Returns:
        Nothing

    
`strip_json_files(source_dir, target_dir, pattern='*.json.gz')`
:   Strip all JSON files of a project by keeping only 'PASS' variants.
    
    JSON files from Illumina's ICA pipeline can be very large because they
    contain any deviation from the reference genome, irrespective of the
    quality of the mutation call. Gzip compressed JSON files with sizes in the
    gigabyte range cannot be processed by JSON packages that read the entire
    file into memory. It is necessary to first reduce the size of JSON files by
    removing all variants that do not meet Illumina's quality criteria.
    
    This function searches `source_dir` recursively for all files matching the
    `file_pattern`. Each of those files is processed and a stripped version
    keeping only variants that PASS Illumina's quality criteria is created. The
    output file has the same name as the input file. The directory structure
    below `source_dir` is replicated in `target_dir`. Output files get the
    suffix '_filtered.json.gz'.‚
    
    Args:
        source_dir: directory where to search for input JSON files
        target_dir: directory where to save the stripped outpout JSON files
        file_pattern: files matching this pattern will be processed
    
    Returns:
        Nothing
    
    Example:
        >>> strip_json_files('../Data/Original', '../Data/Derived')