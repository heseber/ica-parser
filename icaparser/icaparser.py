"""Parser for JSON files from Illumina Connected Annotations pipeline."""

import gzip
import json
import os
import re
import warnings
from collections.abc import Callable
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pkg_resources
from natsort import natsorted
from tqdm.auto import tqdm

_CLINVAR_ORDERED_SIGNIFICANCES = [
    "pathogenic",
    "likely pathogenic",
    "drug response",
    "risk factor",
    "protective",
    "affects",
    "association",
    "uncertain significance",
    "conflicting data from submitters",
    "other",
    "not provided",
    "likely benign",
    "benign",
]

# Adapted from https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl
# Added 'conservative_inframe_deletion':6
# Added 'structural_interaction_variant':1
# Added 'sequence_feature':7
_CONSEQUENCE_PRIORITY = {
    # Added because missing in vcf2maf.pl but contained in Almac's files
    "conservative_inframe_deletion": 6,
    "structural_interaction_variant": 1,
    "sequence_feature": 7,
    # Added because missing in vcr2map.pl but contained in ICA's RNAseq
    # files. In the examples I have seen so far, "transcript_variant" means
    # a SNV, but not always. When investigating individual examples, the specied
    # base change is not in agreement with the reference base of the annotated
    # transcript, so it is unclear whether these changes are reliable. I assign
    # a priority of 6 to such variants and leave it up to the user to keep or
    # remove transcript_variant variants.
    "transcript_variant": 6,
    # from vcf2maf.pl
    "transcript_ablation": 1,  ## A feature ablation whereby the deleted region includes a transcript feature
    "exon_loss_variant": 1,  ## A sequence variant whereby an exon is lost from the transcript
    "splice_donor_variant": 2,  ## A splice variant that changes the 2 base region at the 5' end of an intron
    "splice_acceptor_variant": 2,  ## A splice variant that changes the 2 base region at the 3' end of an intron
    "stop_gained": 3,  ## A sequence variant whereby at least one base of a codon is changed,
    # resulting in a premature stop codon, leading to a shortened transcript
    "frameshift_variant": 3,  ## A sequence variant which causes a disruption of the translational reading
    # frame, because the number of nucleotides inserted or deleted is not a
    # multiple of three
    "stop_lost": 3,  # A sequence variant where at least one base of the terminator codon (stop)
    # is changed, resulting in an elongated transcript
    "start_lost": 4,  ## A codon variant that changes at least one base of the canonical start codon
    "initiator_codon_variant": 4,  ## A codon variant that changes at least one base of the first codon of a
    # transcript
    "disruptive_inframe_insertion": 5,  ## An inframe increase in cds length that inserts one or more codons into the
    # coding sequence within an existing codon
    "disruptive_inframe_deletion": 5,  ## An inframe decrease in cds length that deletes bases from the coding sequence
    # starting within an existing codon
    "inframe_insertion": 5,  ## An inframe non synonymous variant that inserts bases into the coding sequence
    "inframe_deletion": 5,  ## An inframe non synonymous variant that deletes bases from the coding sequence
    "protein_altering_variant": 5,  ## A sequence variant which is predicted to change the protein encoded in the
    # coding sequence
    "missense_variant": 6,  ## A sequence variant, that changes one or more bases, resulting in a different
    # amino acid sequence but where the length is preserved
    "conservative_missense_variant": 6,  ## A sequence variant whereby at least one base of a codon is changed resulting
    # in a codon that encodes for a different but similar amino acid. These
    # variants may or may not be deleterious
    "rare_amino_acid_variant": 6,  ## A sequence variant whereby at least one base of a codon encoding a rare amino
    # acid is changed, resulting in a different encoded amino acid
    "transcript_amplification": 7,  ## A feature amplification of a region containing a transcript
    "splice_region_variant": 8,  ## A sequence variant in which a change has occurred within the region of the
    # splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
    "start_retained_variant": 9,  ## A sequence variant where at least one base in the start codon is changed, but
    # the start remains
    "stop_retained_variant": 9,  ## A sequence variant where at least one base in the terminator codon is changed,
    # but the terminator remains
    "synonymous_variant": 9,  ## A sequence variant where there is no resulting change to the encoded amino
    # acid
    "incomplete_terminal_codon_variant": 10,  ## A sequence variant where at least one base of the final codon of an
    # incompletely annotated transcript is changed
    "coding_sequence_variant": 11,  ## A sequence variant that changes the coding sequence
    "mature_miRNA_variant": 11,  ## A transcript variant located with the sequence of the mature miRNA
    "exon_variant": 11,  ## A sequence variant that changes exon sequence
    "5_prime_UTR_variant": 12,  ## A UTR variant of the 5' UTR
    "5_prime_UTR_premature_start_codon_gain_variant": 12,  # snpEff-specific effect, creating a start codon in 5' UTR
    "3_prime_UTR_variant": 12,  ## A UTR variant of the 3' UTR
    "non_coding_exon_variant": 13,  ## A sequence variant that changes non-coding exon sequence
    "non_coding_transcript_exon_variant": 13,  ## snpEff-specific synonym for non_coding_exon_variant
    "non_coding_transcript_variant": 14,  ## A transcript variant of a non coding RNA gene
    "nc_transcript_variant": 14,  ## A transcript variant of a non coding RNA gene (older alias for
    # non_coding_transcript_variant)
    "intron_variant": 14,  ## A transcript variant occurring within an intron
    "intragenic_variant": 14,  ## A variant that occurs within a gene but falls outside of all transcript
    # features. This occurs when alternate transcripts of a gene do not share
    # overlapping sequence
    "INTRAGENIC": 14,  ## snpEff-specific synonym of intragenic_variant
    "NMD_transcript_variant": 15,  ## A variant in a transcript that is the target of NMD
    "upstream_gene_variant": 16,  ## A sequence variant located 5' of a gene
    "downstream_gene_variant": 16,  ## A sequence variant located 3' of a gene
    "TFBS_ablation": 17,  ## A feature ablation whereby the deleted region includes a transcription factor
    # binding site
    "TFBS_amplification": 17,  ## A feature amplification of a region containing a transcription factor binding
    # site
    "TF_binding_site_variant": 17,  ## A sequence variant located within a transcription factor binding site
    "regulatory_region_ablation": 17,  ## A feature ablation whereby the deleted region includes a regulatory region
    "regulatory_region_amplification": 17,  ## A feature amplification of a region containing a regulatory region
    "regulatory_region_variant": 17,  ## A sequence variant located within a regulatory region
    "regulatory_region": 17,  ## snpEff-specific effect that should really be regulatory_region_variant
    "feature_elongation": 18,  ## A sequence variant that causes the extension of a genomic feature, with
    # regard to the reference sequence
    "feature_truncation": 18,  ## A sequence variant that causes the reduction of a genomic feature, with
    # regard to the reference sequence
    "intergenic_variant": 19,  ## A sequence variant located in the intergenic region, between genes
    "intergenic_region": 19,  ## snpEff-specific effect that should really be intergenic_variant
    "": 20,
}

_BIOTYPE_PRIORITY = {
    "prime3_overlapping_ncrna": 5,  ## added by Henrik, probably a typo in Almac's annotation,
    # should be 3prime_overlapping_ncrna which has prio 5
    "protein_coding": 1,  ## Contains an open reading frame (ORF)
    "mRNA": 1,  ## Newer versions of ICA use "mRNA" instead of "protein_coding"
    "LRG_gene": 2,  ## Gene in a "Locus Reference Genomic" region known to have disease-related
    # sequence variations
    "IG_C_gene": 2,  ## Immunoglobulin (Ig) variable chain genes imported or annotated according
    # to the IMGT
    "IG_D_gene": 2,  ## Immunoglobulin (Ig) variable chain genes imported or annotated according
    #  to the IMGT
    "IG_J_gene": 2,  ## Immunoglobulin (Ig) variable chain genes imported or annotated according
    # to the IMGT
    "IG_LV_gene": 2,  ## Immunoglobulin (Ig) variable chain genes imported or annotated according
    # to the IMGT
    "IG_V_gene": 2,  ## Immunoglobulin (Ig) variable chain genes imported or annotated according
    # to the IMGT
    "TR_C_gene": 2,  ## T-cell receptor (TcR) genes imported or annotated according to the IMGT
    "TR_D_gene": 2,  ## T-cell receptor (TcR) genes imported or annotated according to the IMGT
    "TR_J_gene": 2,  ## T-cell receptor (TcR) genes imported or annotated according to the IMGT
    "TR_V_gene": 2,  ## T-cell receptor (TcR) genes imported or annotated according to the IMGT
    "miRNA": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "snRNA": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "snoRNA": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "ribozyme": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "tRNA": 3,  ##Added by Y. Boursin
    "sRNA": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "scaRNA": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "rRNA": 3,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "scRNA": 3,  ## Non-coding RNA predicted using sequences from Rfam and miRBase
    "lincRNA": 3,  ## Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily
    # conserved, intergenic regions
    "lncRNA": 3,  ## Replaces 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA,
    # lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic
    # and sense_overlapping
    "bidirectional_promoter_lncrna": 3,  ## A non-coding locus that originates from within the promoter region of a
    # protein-coding gene, with transcription proceeding in the opposite
    # direction on the other strand
    "bidirectional_promoter_lncRNA": 3,  ## A non-coding locus that originates from within the promoter region of a
    # protein-coding gene, with transcription proceeding in the opposite
    # direction on the other strand
    "known_ncrna": 4,
    "vaultRNA": 4,  ## Short non coding RNA genes that form part of the vault ribonucleoprotein
    # complex
    "macro_lncRNA": 4,  ## unspliced lncRNAs that are several kb in size
    "Mt_tRNA": 4,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "Mt_rRNA": 4,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "antisense": 5,  ## Has transcripts that overlap the genomic span (i.e. exon or introns) of a
    # protein-coding locus on the opposite strand
    "antisense_RNA": 5,  ## Alias for antisense (Y. Boursin)
    "sense_intronic": 5,  ## Long non-coding transcript in introns of a coding gene that does not overlap
    # any exons
    "sense_overlapping": 5,  ## Long non-coding transcript that contains a coding gene in its intron on the
    # same strand
    "3prime_overlapping_ncrna": 5,  ## Transcripts where ditag and/or published experimental data strongly supports
    # the existence of short non-coding transcripts transcribed from the 3'UTR
    "3prime_overlapping_ncRNA": 5,  ## Transcripts where ditag and/or published experimental data strongly supports
    # the existence of short non-coding transcripts transcribed from the 3'UTR
    "misc_RNA": 5,  ## Non-coding RNA predicted using sequences from RFAM and miRBase
    "non_coding": 5,  ## Transcript which is known from the literature to not be protein coding
    "regulatory_region": 6,  ## A region of sequence that is involved in the control of a biological process
    "disrupted_domain": 6,  ## Otherwise viable coding region omitted from this alternatively spliced
    # transcript because the splice variation affects a region coding for a
    # protein domain
    "processed_transcript": 6,  ## Doesn't contain an ORF
    "TEC": 6,  ## To be Experimentally Confirmed. This is used for non-spliced EST clusters
    # that have polyA features. This category has been specifically created for
    # the ENCODE project to highlight regions that could indicate the presence of
    # protein coding genes that require experimental validation, either by 5'
    # RACE or RT-PCR to extend the transcripts, or by confirming expression of
    # the putatively-encoded peptide with specific antibodies
    "TF_binding_site": 7,  ## A region of a nucleotide molecule that binds a Transcription Factor or
    # Transcription Factor complex
    "CTCF_binding_site": 7,  ## A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG,
    # bound by CCCTF-binding factor
    "promoter_flanking_region": 7,  ## A region immediately adjacent to a promoter which may or may not contain
    # transcription factor binding sites
    "enhancer": 7,  ## A cis-acting sequence that increases the utilization of (some) eukaryotic
    # promoters, and can function in either orientation and in any location
    # (upstream or downstream) relative to the promoter
    "promoter": 7,  ## A regulatory_region composed of the TSS(s) and binding sites for
    # TF_complexes of the basal transcription machinery
    "open_chromatin_region": 7,  ## A DNA sequence that in the normal state of the chromosome corresponds to an
    # unfolded, un-complexed stretch of double-stranded DNA
    "retained_intron": 7,  ## Alternatively spliced transcript believed to contain intronic sequence
    # relative to other, coding, variants
    "nonsense_mediated_decay": 7,  ## If the coding sequence (following the appropriate reference) of a transcript
    # finishes >50bp from a downstream splice site then it is tagged as NMD. If
    # the variant does not cover the full reference coding sequence then it is
    # annotated as NMD if NMD is unavoidable i.e. no matter what the exon
    # structure of the missing portion is the transcript will be subject to NMD
    "non_stop_decay": 7,  ## Transcripts that have polyA features (including signal) without a prior stop
    # codon in the CDS, i.e. a non-genomic polyA tail attached directly to the
    # CDS without 3' UTR. These transcripts are subject to degradation
    "ambiguous_orf": 7,  ## Transcript believed to be protein coding, but with more than one possible
    # open reading frame
    "pseudogene": 8,  ## Have homology to proteins but generally suffer from a disrupted coding
    # sequence and an active homologous gene can be found at another locus.
    # Sometimes these entries have an intact coding sequence or an open but
    # truncated ORF, in which case there is other evidence used (for example
    # genomic polyA stretches at the 3' end) to classify them as a pseudogene.
    # Can be further classified as one of the following
    "processed_pseudogene": 8,  ## Pseudogene that lack introns and is thought to arise from reverse
    # transcription of mRNA followed by reinsertion of DNA into the genome
    "polymorphic_pseudogene": 8,  ## Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains
    # the gene is translated
    "retrotransposed": 8,  ## Pseudogene owing to a reverse transcribed and re-inserted sequence
    "translated_processed_pseudogene": 8,  ## Pseudogenes that have mass spec data suggesting that they are also translated
    "translated_unprocessed_pseudogene": 8,  ## Pseudogenes that have mass spec data suggesting that they are also translated
    "transcribed_processed_pseudogene": 8,  ## Pseudogene where protein homology or genomic structure indicates a pseudogene,
    # but the presence of locus-specific transcripts indicates expression
    "transcribed_unprocessed_pseudogene": 8,  ## Pseudogene where protein homology or genomic structure indicates a pseudogene,
    # but the presence of locus-specific transcripts indicates expression
    "transcribed_unitary_pseudogene": 8,  ##Pseudogene where protein homology or genomic structure indicates a pseudogene,
    # but the presence of locus-specific transcripts indicates expression
    "unitary_pseudogene": 8,  ## A species specific unprocessed pseudogene without a parent gene, as it has an
    # active orthologue in another species
    "unprocessed_pseudogene": 8,  ## Pseudogene that can contain introns since produced by gene duplication
    "Mt_tRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "tRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "snoRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "snRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "scRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "rRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "misc_RNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "miRNA_pseudogene": 8,  ## Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
    "IG_C_pseudogene": 8,  ## Inactivated immunoglobulin gene
    "IG_D_pseudogene": 8,  ## Inactivated immunoglobulin gene
    "IG_J_pseudogene": 8,  ## Inactivated immunoglobulin gene
    "IG_V_pseudogene": 8,  ## Inactivated immunoglobulin gene
    "TR_J_pseudogene": 8,  ## Inactivated immunoglobulin gene
    "TR_V_pseudogene": 8,  ## Inactivated immunoglobulin gene
    "artifact": 9,  ## Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
    "": 10,
}

# Read, sort and rank the VEP consequences from a file in the data directory
vep_csq = pd.read_csv(
    pkg_resources.resource_stream(__name__, "data/vep_consequences.txt"), sep="\t"
)
vep_csq = vep_csq.set_index("SO term")
vep_csq["rank"] = list(range(1, 1 + vep_csq.shape[0]))
vep_csq["priority"] = [_CONSEQUENCE_PRIORITY[x] for x in vep_csq.index]
vep_csq = vep_csq.sort_values(["priority", "rank"])

# Read the oncogene / tumor suppressor gene definition from a file in the data directory
gene_type = pd.read_csv(
    pkg_resources.resource_stream(
        __name__, "data/cancer_genes_union_TCGA_COSMIC_manual.tsv"
    ),
    sep="\t",
    index_col="gene_symbol",
)


def strip_json_file(ifname: str, ofname: str) -> None:
    """Reduce the JSON file size by keeping only 'PASS' variants.

    JSON files from Illumina's ICA pipeline can be very large because they
    contain any deviation from the reference genome, irrespective of the
    quality of the mutation call. Gzip compressed JSON files with sizes in the
    gigabyte range cannot be processed by JSON packages that read the entire
    file into memory. It is necessary to first reduce the size of JSON files by
    removing all variants that do not meet Illumina's quality criteria.

    This function reads a single JSON file and creates a single JSON outpout
    file by removing all variants that do not pass Illumina's quality criteria.

    Args:
        ifname: name of the input file.
        ofname: name of the output file.

    Returns:
        None.
    """
    if ifname == ofname:
        raise ValueError("ifname and ofname must be different")
    odname = os.path.dirname(ofname)
    os.makedirs(odname, exist_ok=True)
    is_header_line = True
    is_position_line = False
    is_first_position_line = False
    is_gene_line = False
    is_first_gene_line = False
    gene_section_line = '],"genes":['
    end_line = "]}"
    with gzip.open(ifname, "rt") as fin:
        with gzip.open(ofname, "wt") as fout:
            for line in fin:
                trim_line = line.strip()
                if is_header_line:
                    print(trim_line, file=fout)
                    is_header_line = False
                    is_position_line = True
                    is_first_position_line = True
                    continue
                if trim_line == gene_section_line:
                    print("", trim_line, sep="\n", file=fout)
                    is_gene_line = True
                    is_first_gene_line = True
                    is_position_line = False
                    continue
                elif trim_line == end_line:
                    print("", trim_line, sep="\n", file=fout)
                    break
                else:
                    if is_position_line:
                        ## remove the trailing ',' if there is
                        position_json = trim_line.rstrip(",")
                        position = json.loads(position_json)
                        if "PASS" in position["filters"]:
                            if is_first_position_line:
                                finish_prev_line = ""
                            else:
                                finish_prev_line = ",\n"
                            print(
                                finish_prev_line,
                                position_json,
                                sep="",
                                end="",
                                file=fout,
                            )
                            is_first_position_line = False
                    if is_gene_line:
                        ## remove the trailing ',' if there is
                        gene_json = trim_line.rstrip(",")
                        if is_first_gene_line:
                            finish_prev_line = ""
                        else:
                            finish_prev_line = ",\n"
                        print(finish_prev_line, gene_json, sep="", end="", file=fout)
                        is_first_gene_line = False


def strip_json_files(
    source_dir: str, target_dir: str, pattern: str = "*.json.gz"
) -> None:
    """Strip all JSON files of a project by keeping only 'PASS' variants.

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
    suffix '_filtered.json.gz'.

    Args:
        source_dir: directory where to search for input JSON files.
        target_dir: directory where to save the stripped outpout JSON files.
        pattern: files matching this pattern will be processed.

    Returns:
        None.

    Examples:
        >>> strip_json_files('../Data/Original', '../Data/Derived')
    """
    source_dir = Path(source_dir).absolute()
    target_dir = Path(target_dir).absolute()
    if source_dir == target_dir:
        raise ValueError("source and target directory must be different")
    for path in tqdm(list(source_dir.rglob(pattern))):
        ifname = path.as_posix()
        ofname = ifname.replace(source_dir.as_posix(), target_dir.as_posix(), 1)
        ofname = ofname.removesuffix(".json.gz") + "_filtered.json.gz"
        tqdm.write("Stripping " + ifname)
        strip_json_file(ifname, ofname)


def get_dna_json_files(
    base_dir: str, pattern: str = "*MergedVariants_Annotated_filtered.json.gz"
) -> list:
    """Find DNA annotation JSON files in or below `base_dir`.

    Searches for ICA DNA annotation JSON files in and below `base_dir`.
    All file names matching `pattern` are returned.

    Args:
        base_dir: base directory of directory subtree where to search
                  for DNA annotation JSON files.
        pattern:  files names matching this pattern are returned.

    Returns:
        file names.
    """
    files = [x.as_posix() for x in Path(base_dir).rglob(pattern)]
    files = sorted(files)
    return files


def get_header(file: str) -> dict:
    """Extract the header element from a ICA JSON file.

    Args:
        file: name of the ICA JSON file.

    Returns:
        header from the JSON file.
    """
    with gzip.open(file, "rt") as f:
        line = next(f)
    trim_line = line.strip()
    header = trim_line[10:-14]
    header = json.loads(header)
    return header


def get_sample(file: str, suffix: str = "(-D[^.]*)?\\.bam") -> str:
    """Extract the sample name from a ICA JSON file.

    Args:
        file: name of the ICA JSON file.
        suffix: regular expression to remove from the sample name in the JSON
            file. Defaults to '(-D[^.]*)?\\.bam'.

    Returns:
        name of the sample annotated in the JSON file.
    """
    header = get_header(file)
    sample = header["samples"][0]
    sample = re.sub(suffix, "", sample)
    return sample


def get_header_scalars(file: str) -> pd.DataFrame:
    """Extract a table with all scalar attributes from the JSON header.

    Args:
        file: name of the ICA JSON file.

    Returns:
        table of scalar attributes and their values.
    """
    header = get_header(file)
    sample = get_sample(file)
    df = pd.DataFrame([x for x in header.items() if not isinstance(x[1], list)])
    df = df.rename(columns={0: "Attribute", 1: sample})
    return df


def get_data_sources(file: str) -> pd.DataFrame:
    """Extract a table with annotation data sources from the JSON header.

    Args:
        file: name of the ICA JSON file.

    Returns:
        table with annotation data sources and their versions.
    """
    header = get_header(file)
    df = pd.DataFrame(header["dataSources"])
    return df


def get_pipeline_metadata(files: list) -> pd.DataFrame:
    """Extract a table with metadata annotation pipeline run from the JSON header.

    Args:
        files: names of the ICA JSON files.

    Returns:
        table with metadata of pipeline runs.
    """
    df = get_header_scalars(files[0])
    for fname in files[1:]:
        df = pd.merge(df, get_header_scalars(fname), on="Attribute", how="outer")
    df = df.set_index("Attribute").transpose()
    df.index.rename("Sample", inplace=True)
    return df


def _get_branch(file: str, key: str) -> list:
    """Extract a subsection (branch) of the ICA JSON file.

    This is a helper function which is not supposed to be used by the end user.

    Args:
        file: name of the ICA JSON file.
        key:  name of the top level element to extract (`positions` or `genes`).

    Returns:
        list of positions or genes.

    """
    with gzip.open(file, "rt") as f:
        data = json.load(f)
    if key in data.keys():
        return data[key]
    else:
        return list()


def cleanup_cosmic(positions: list) -> list:
    """Remove Cosmic entries with alleles not matching the variant alleles.

    ICA attaches Cosmic entries to variants based on position only, which
    leads to wrong assignments of Cosmic entries to variants. This function
    removes all Cosmic entries from a variant for which reference and altered
    alleles do not match those of the variant.

    Filtering is done in place.

    Args:
        positions: list of positions to clean up.

    Returns:
        list of positions with cleaned up Cosmic entries.
    """
    for p in positions:
        for v in p["variants"]:
            if "cosmic" not in v:
                continue
            cosmic = []
            for c in v.get["cosmic"]:
                # if c['refAllele'] == variant['refAllele'] and c['altAllele'] == variant['altAllele']:
                if c.get("isAlleleSpecific", False):
                    cosmic.append(c)
            if len(cosmic) > 0:
                v["cosmic"] = cosmic
            else:
                del v["cosmic"]
    return positions


def get_positions(
    file: str, variant_filters: list = [], transcript_filters: list = []
) -> list:
    """Extract all positions from a ICA JSON file.

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
        filtered positions from file.

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
    """
    positions = _get_branch(file, "positions")
    positions = add_gene_types(positions)
    for vf in variant_filters:
        positions = filter_positions_by_variants(positions, vf)
    for tf in transcript_filters:
        positions = filter_positions_by_transcripts(positions, tf)
    sample_id = get_sample(file)
    for p in positions:
        p["samples"][0]["sampleId"] = sample_id
    return positions


def get_multi_sample_positions(files: list, *args: object, **kwargs: object) -> list:
    """Extract all positions for a set of ICA JSON files.

    The sample id is stored as an additional new attribute of the
    `samples` element of a position. The `samples` element is a list,
    although ICA usually only creates single sample JSON files.

    Args:
        files: names of the ICA JSON files.
        args: extra arguments forwarded to get_positions().
        kwargs: extra named arguments forwarded to get_positions().

    Returns:
        filtered positions from all files.

    Examples:
        >>> import icaparser as icap
        >>> positions = icap.get_multi_sample_positions(json_files)
        >>> print(positions[0]['samples'][0]['sampleId'])
    """
    positions = []
    for f in tqdm(files):
        positions += get_positions(f, *args, **kwargs)
    return positions


def get_genes(file: str) -> list:
    """Extract gene annotation from a ICA JSON file.

    The `genes` section of ICA JSON files is optional. If this section
    is not included in the file, an empty list is returned.

    Args:
        file: name of the ICA JSON file.

    Returns:
        gene annotations.
    """
    return _get_branch(file, "genes")


def get_position_by_coordinates(
    positions: list, chromosome: str, position: int
) -> dict:
    """Extract a particular position from a position list.

    Args:
        positions: list of input positions.
        chromosome: name of the chromosome.
        position: numeric position on the chromosome.

    Returns:
        the position for the specified chromosome and numeric position.

    Examples:
        >>> import icaparser as icap
        >>> icap.get_position_by_coordinates(positions, 'chr1', 204399064)
    """

    def filter_func(x):
        return x["chromosome"] == chromosome and x["position"] == int(position)

    return list(filter(filter_func, positions))[0]


def get_max_af(variant: dict, source: str, cohorts: list = None) -> float:
    """Get the maximum allele frequency for a particular annotation source.

    Get the maximum allele frequency across all cohorts annotated by the
    annotation source.

    Args:
        variant: the variant to investigate.
        source: the annotation source to use, for example 'gnomad' or
            'gnomadExome' or 'oneKg'.
        cohorts: subpopulations to include; include all if None.

    Returns:
        the maximum allele frequency.

    Examples:
        >>> import icaparser as icap
        >>> icap.get_max_af(variant, 'gnomad')
    """
    if source not in variant:
        return 0
    max_af = 0
    for k, v in variant[source].items():
        if k.endswith("Af"):
            if not cohorts or k.replace("Af", "") in cohorts:
                max_af = max(max_af, v)
    return max_af


def get_gnomad_max_af(
    variant: dict, cohorts: list = ["afr", "amr", "eas", "nfe", "sas"]
) -> float:
    """Get the maximum allele frequency for gnomAD.

    Get the maximum allele frequences across all major cohorts annotated
    by gnomAD, excluding bottleneck populations (Ashkenazy Jews and Finish)
    and _other_.

    Args:
        variant: the variant to investigate.
        cohorts: subpopulations to include.

    Returns:
        maximum GnomAD allele frequency.
    """
    return get_max_af(variant, "gnomad", cohorts)


def get_gnomad_exome_max_af(
    variant: dict, cohorts: list = ["afr", "amr", "eas", "nfe", "sas"]
) -> float:
    """Get the maximum allele frequency for gnomAD Exome.

    Get the maximum allele frequences across all major cohorts annotated
    by gnomAD, Exome  excluding bottleneck populations (Ashkenazy Jews and
    Finish) and _other_.

    Args:
        variant: the variant to investigate.
        cohorts: subpopulations to include.

    Returns:
        maximum GnomAD Exome allele frequency.
    """
    return get_max_af(variant, "gnomadExome", cohorts)


def get_onekg_max_af(variant: dict) -> float:
    """Get the maximum allele frequency for the 1000 Genomes Project.

    Get the maximum allele frequences across all cohorts annotated
    by the 1000 Genomes Project.

    Args:
        variant: the variant to investigate.

    Returns:
        maximum 1000 genomes allele frequency.
    """
    return get_max_af(variant, "oneKg")


def get_cosmic_max_sample_count(
    variant: dict, only_allele_specific: bool = True
) -> int:
    """Get the maximum sample count for all Cosmic annotations of a variant.

    A variant can have no, one or multiple associated Cosmic identifiers. This
    function returns the maximum sample count of all Cosmic identifiers. For
    each Cosmic identifier, sample numbers are summed up across all indications.
    Returns 0 if no Cosmic identifier exists for this variant.

    The 'only_allele_specific' argument is used to exclude Cosmic entries that
    annotate the same chromosomal location but an allele that is different from
    the allele of the annotated variant. ICA annotates a variant with all Cosmic
    entries for that chromosomal location, irrespective of alleles. When
    counting Cosmic samples, this leads to an overestimation of Cosmic sample
    counts for a particular variant. Therefore, 'only_allele_specific' is True
    by default to count only samples from Cosmic entries with matching alleles.
    Occasionally, it may be desired, though, to count all samples with mutations
    at a given position, irrespective of allele. For example, several different
    alleles at a functional site of a gene can lead to function-disrupting
    mutations, so we want to get the maximum sample count for any allele at that
    position. One might also think of adding the sample counts for all Cosmic
    entries annotating a variant, but this does not work  due to redundancy of
    Cosmic entries. Older Cosmic versions often included the same sample in
    different Cosmic entries. And newer Cosmic versions often have multiple
    entries for an allele, one for each transcript variant, with the same
    underlying samples.

    Args:
        variant:              the variant to investigate
        only_allele_specific: consider only cosmic entries with alleles
                              matching the allele of the annotated variant

    Returns:
        maximum cosmic sample count
    """
    max_count = 0
    for cosmic_entry in variant.get("cosmic", []):
        if only_allele_specific and not cosmic_entry.get("isAlleleSpecific", False):
            continue
        max_count = max(max_count, cosmic_entry.get("sampleCount", 0))
    return max_count


def get_clinvar_max_significance(
    variant: dict, ordered_significances: list = _CLINVAR_ORDERED_SIGNIFICANCES
) -> str:
    """Get the maximum signifinance for all ClinVar annotations of a variant.

    Args:
        variant: the variant to investigate.
        ordered_significances: ranked order of ClinVar significances.

    Returns:
        ClinVar significance of highest rank for the variant.
    """
    if "clinvar" not in variant:
        return "none"
    significance_ranks = dict(
        zip(ordered_significances, list(range(len(ordered_significances))))
    )
    best_rank = max(significance_ranks.values()) + 1
    max_significance = "none"
    for clinvar in variant["clinvar"]:
        for significance in clinvar["significance"]:
            if significance not in significance_ranks:
                continue
            rank = significance_ranks[significance]
            if rank < best_rank:
                best_rank = rank
                max_significance = significance
    return max_significance


def common_variant_filter(variant: dict, max_af: float = 0.001) -> bool:
    """Get a variant filter based on GnomAD, GnomAd Exome, and 1000 Genomes.

    Returns True if none of the maximum allele frequencies from GnomAD, GnomAD
    exomes and 1000 genomes is greater than `max_af`. The default value of 0.1 %
    for the maximum allele frequency corresponds to that of the AACR GENIE
    project.

    Args:
        variant: the variant to investigate.
        max_af: the maximum allele frequency threshold.

    Returns:
        True if this is not a common variant.
    """
    passed = (
        get_gnomad_max_af(variant) <= max_af
        and get_gnomad_exome_max_af(variant) <= max_af
        and get_onekg_max_af(variant) <= max_af
    )
    return passed


def _filter_items(
    items: list, filter_func: Callable[[dict], bool], sub_items_name: str
) -> list:
    """Filter items based on a filter for subitems.

    Args:
        items: list of items to filter.
        filter_func: function returning a bool for each sub item.
        sub_items_name: each item is a dict, and this is the key
            in this item dict for which the values are a list of sub items

    Returns:
        filtered items.
    """
    filtered_items = []
    for item in items:
        if sub_items_name not in item:
            continue
        sub_items = item[sub_items_name]
        filtered_sub_items = list(filter(filter_func, sub_items))
        if filtered_sub_items:
            filtered_item = deepcopy(item)
            filtered_item[sub_items_name] = filtered_sub_items
            filtered_items.append(filtered_item)
    return filtered_items


def filter_positions_by_variants(
    positions: list, filter_func: Callable[[dict], bool]
) -> list:
    """Filter positions based on a filter function for variants.

    Apply a filter function to all variants of each position.
    Variants not passing the filter are removed from a position.
    Positions without any variants passing the filter are removed
    from the returned list.

    Args:
        positions: list of positions to filter.
        filter_func: function taking a variant and returning a bool.
            True means to keep the variant.

    Returns:
        filtered positions.

    Examples:
        >>> import icaparser as icap
        >>> max_af = 0.001
        >>> is_not_common_variant = lambda x: icap.common_variant_filter(x, max_af)
        >>> non_common_positions = icap.filter_positions_by_variants(
                positions,
                is_not_common_variant
            )
    """
    return _filter_items(positions, filter_func, "variants")


def filter_variants_by_transcripts(
    variants: list, filter_func: Callable[[dict], bool]
) -> list:
    """Filter variants based on a filter function for transcripts.

    Apply a filter function to all transcripts of each variant.
    Transcripts not passing the filter are removed from a variant.
    Variants without any transcripts passing the filter are removed
    from the returned list.

    Args:
        variants: list of variants to filter.
        filter_func: function taking a transcript and returning a bool.
            True means to keep the transcript.

    Returns:
        filtered variants.
    """
    return _filter_items(variants, filter_func, "transcripts")


def filter_positions_by_transcripts(
    positions: list, filter_func: Callable[[dict], bool]
) -> list:
    """Filter positions based on a filter function for transcripts.

    Apply a filter function to all transcripts of each position. Transcripts not
    passing the filter are removed from the variants of a position. Variants
    without any transcript left are removed from a position. Positions without
    any variants left are removed from the returned list of positions.

    Args:
        positions: list of positions to filter.
        filter_func: function taking a transcript and returning a bool.
            True means to keep the transcript.

    Returns:
        filtered positions.

    Examples:
        >>> is_canonical_transcript = lambda x: x.get('isCanonical', False)
        >>> canonical_positions = icap.filter_positions_by_transcripts(
                non_common_positions,
                is_canonical_transcript
            )
    """
    filtered_positions = []
    for position in positions:
        variants = position["variants"]
        filtered_variants = filter_variants_by_transcripts(variants, filter_func)
        if filtered_variants:
            filtered_position = deepcopy(position)
            filtered_position["variants"] = filtered_variants
            filtered_positions.append(filtered_position)
    return filtered_positions


def get_clinvar(variant: dict) -> pd.DataFrame:
    """Get a table of all ClinVar annotations for a variant.

    Args:
        variant: the variant to investigate.

    Returns:
        table with ClinVar annotations.
    """
    if "clinvar" not in variant:
        return pd.DataFrame()
    return pd.DataFrame(variant["clinvar"])


def get_biotype_priority(biotype: str) -> int:
    """Get the numeric priority of a biotype.

    The numeric priority of a biotype that is returned by this function is the
    same as defined by [vcf2maf.pl by MSKCC](
    https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl ). Biotypes are
    'protein_coding', 'LRG_gene', ,'miRNA', ...

    Args:
        biotype: the biotype for which the priority is to be returned.

    Returns:
        the priority, smaller values mean higher priority.

    """
    lowest_prio = max(_BIOTYPE_PRIORITY.values())
    if biotype not in _BIOTYPE_PRIORITY.keys():
        warnings.warn(
            '"'
            + biotype
            + '" is not a known biotype and will be set to the lowest priority of '
            + str(lowest_prio)
            + "."
        )
    return _BIOTYPE_PRIORITY.get(biotype, lowest_prio)


def get_consequences(transcript: dict) -> list:
    """Get a list of consequences for a transcript.

    A list of consequences of a variant for a transcript is returned. If any of
    the annotated consequences is a combination of single consequences,
    separated by ampersands (&) or commas, the consequence is split into
    single consequences.

    Args:
        transcript: the transcript for which the consequences are to be returned.

    Returns:
        the consequences, a list of strings.
    """
    consequences = transcript.get("consequence", [])
    consequences = [x for c in consequences for x in re.split("&,", c)]
    return consequences


def get_vep_rank_for_consequence(consequence: str) -> int:
    """Get the numeric rank of a VEP consequence term.

    The numeric rank of a consequence is the position of the consequence in this
    [list of consequences for the Variant Effect Predictor VEP](
    https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    ).

    Args:
        consequence: the consequence term of the variant.

    Returns:
        the rank of the consequence, smaller values mean higher rank.
    """
    if consequence not in vep_csq.index:
        return vep_csq["rank"].max() + 1
    return vep_csq.loc[consequence, "rank"]


def get_vep_priority_for_consequence(consequence: str):
    """Get the numeric priority of a VEP consequence term.

      The numeric priority of a consequence that is returned by this function is
      the same as defined by [vcf2maf.pl of MSKCC](
      https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl ).

    Args:
          consequence: the consequence term of the variant.

      Returns:
          the priority of the consequence, smaller values mean higher priority.
    """
    if consequence not in vep_csq.index:
        return vep_csq["priority"].max() + 1
    return vep_csq.loc[consequence, "priority"]


def get_vep_consequence_for_rank(rank: int) -> str:
    """Get the VEP consequence term of a numeric rank.

    Args:
        rank: the numeric rank of the consequence term.

    Returns:
        the consequence.
    """
    if rank > vep_csq.shape[0]:
        return "unknown"
    else:
        return list(vep_csq[vep_csq["rank"] == rank].index)[0]


def get_strongest_vep_consequence_rank(transcript: dict) -> int:
    """Get the strongest rank of VEP consequences for a transcript.

    Get the strongest numeric rank of all VEP consequences for a transcript.
    Smaller ranks mean stronger impact.

    The priority of consequences is taken into account first. So if two
    consequences have different priorities, the consequence with the higher
    priority (lower priority number) will be used, and the rank for this
    consequence will be returned. If there are multiple consequences with the
    same priority, the lowest (strongest) rank will be returned.

    For clarification: ranks are unique, i.e. all VEP consequences ordered as
    listed on the VEP documentation page get the row number of this table
    assigned as rank.

    However, several consequences can have the same priority (e.g., stop gained
    and frameshift have the same priority). Priorities are copied from
    [vcf2maf.pl of MSKCC](
    https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl ).

    Args:
        transcript: the transcript to investigate.

    Returns:
        the rank of the VEP consequence with strongest impact.
    """
    if "consequence" not in transcript:
        return vep_csq["rank"].max() + 1
    consequences = get_consequences(transcript)
    df = vep_csq.loc[consequences]
    strongest_prio = df.priority.min()
    min_rank = df[df.priority == strongest_prio]["rank"].min()
    return min_rank


def get_strongest_vep_consequence_priority(transcript: dict) -> int:
    """Get the strongest priority of VEP consequence for a transcript.

    Get the strongest numeric priority of all VEP consequences
    for a transcript. Smaller numeric priorities mean stronger impact.

    Args:
        transcript: the transcript to investigate.

    Returns:
        the strongest numeric priority, smaller values mean higher priority.
    """
    if "consequence" not in transcript:
        return vep_csq["priority"].max() + 1
    consequences = get_consequences(transcript)
    min_prio = vep_csq.loc[consequences, "priority"].min()
    return min_prio


def get_strongest_vep_consequence_name(transcript: dict) -> str:
    """Get the name of the strongest VEP consequence for a transcript.

    Args:
        transcript: the transcript to investigate.

    Returns:
        the consequence.
    """
    rank = get_strongest_vep_consequence_rank(transcript)
    return get_vep_consequence_for_rank(rank)


def get_gene_type(gene_symbol: str) -> str:
    """Get the gene type (oncogene, tsg, mixed) for a gene.

    Args:
        gene_symbol: the gene symbol of the gene.

    Returns:
        the gene type.
    """
    if gene_symbol in gene_type.index:
        return gene_type.loc[gene_symbol, "gene_type"]
    else:
        return ""


def get_mutation_table_for_position(position: dict) -> pd.DataFrame:
    """Get an annotated table of all transcripts for a single position.

    Returns an annotated table of all transcripts that are affected
    by a mutation at a position.

    Args:
        position: the position to investigate.

    Returns:
        table of annotated mutations and affected transcripts.
    """
    vep_csq_categories = list(vep_csq.index)
    vep_csq_categories.reverse()
    rows = []
    for variant_idx, variant in enumerate(position["variants"]):
        # Get the variant frequency for this variant
        variant_frequencies = position["samples"][0].get("variantFrequencies", [])
        if variant_idx < len(variant_frequencies):
            variant_frequency = variant_frequencies[variant_idx]
        else:
            variant_frequency = np.nan
        # Create the table rows
        for transcript in variant["transcripts"]:
            row = {
                "sample": position["samples"][0]["sampleId"],
                "chromosome": position["chromosome"],
                "position": position["position"],
                "genotype": position.get("samples")[0].get("genotype"),
                "variantFrequency": variant_frequency,
                "hgnc": transcript.get("hgnc", ""),
                "source": transcript.get("source", ""),
                "geneId": transcript.get("geneId", ""),
                "transcriptId": transcript.get("transcript", ""),
                "proteinId": transcript.get("proteinId", ""),
                "hgvsc": transcript.get("hgvsc", ""),
                "hgvsp": transcript.get("hgvsp", ""),
                "isCanonical": transcript.get("isCanonical", False),
                "vid": variant["vid"],
                "hgvsg": variant.get("hgvsg", pd.NA),
                "begin": variant["begin"],
                "end": variant["end"],
                "refAllele": variant["refAllele"],
                "altAllele": variant["altAllele"],
                "variantType": variant["variantType"],
                "bioType": transcript.get("bioType", ""),
                "geneType": transcript.get("geneType", ""),
                "mutationStatus": transcript.get("mutationStatus", ""),
                "codons": transcript.get("codons", ""),
                "aminoAcids": transcript.get("aminoAcids", ""),
                "cdnaPos": transcript.get("cdnaPos", ""),
                "cdsPos": transcript.get("cdsPos", ""),
                "exons": transcript.get("exons", ""),
                "proteinPos": transcript.get("proteinPos", ""),
                "consequence": transcript.get("consequence", []),
                "cosmicSampleCount": get_cosmic_max_sample_count(variant),
                "maxGnomadAf": get_gnomad_max_af(variant),
                "maxGnomadExomeAf": get_gnomad_exome_max_af(variant),
                "maxOneKG": get_onekg_max_af(variant),
            }
            rows.append(row)
    df = pd.DataFrame(rows)
    return df


def get_mutation_table_for_positions(
    positions: list, hide_progress: bool = False
) -> pd.DataFrame:
    """Get an annotated table of all transcripts for all positions.

    Returns an annotated table of all transcripts that are affected
    by a mutation at any of the positions.

    Args:
        positions: the positions to investigate.

    Returns:
        table of annotated mutations and affected transcripts.
    """

    if len(positions) == 0:
        return pd.DataFrame()
    mutation_tables = []
    for p in tqdm(positions, disable=hide_progress, desc="Positions"):
        mutation_tables.append(get_mutation_table_for_position(p))
    if len(mutation_tables) == 0:
        return pd.DataFrame()
    mutation_table = pd.concat(mutation_tables, ignore_index=True)
    return mutation_table


def _get_mutation_table_for_single_file(
    file: str,
    max_af: float = 0.001,
    min_vep_consequence_priority: int = 6,
    min_cosmic_sample_count: int = 0,
    only_canonical: bool = False,
    extra_variant_filters: list = [],
    extra_transcript_filters: list = [],
) -> pd.DataFrame:
    """Get an annotated table of all filtered transcripts from a ICA JSON file.

    Load all positions from a single ICA JSON file and filter them. Positions
    having any remaining variants and transcripts passing the filter are
    returned as an annotated table.

    Args:
        file:   an ICA JSON file
        max_af: maximum allele frequency for gnomAD, gnomAD Exome and 1000 Genomes.
                Only variants with maximum allele frequencies below this threshold
                will be returned.
        min_vep_consequence_priority: only transcripts with a minimum VEP
                consequence priority not larger than this threshold will be retained.
                Consequences with priorities <= 6 change the protein sequence,
                consequences with priorities > 6 do not change the protein
                sequence.
        min_cosmic_sample_count: only variants with a maximum cosmic sample count
                not lower than this threshold will be retained
        only_canonical: if true, only canonical transcripts will be retained
        extra_variant_filters: any additional filters to apply to variants.
                Filters shall return True to keep a variant.
        extra_transcript_filters: any additional filters to apply to transcripts.
                Filters shall return True to keep a transcript.

    Returns:
        table of annotated mutations and affected transcripts.

    Examples:
        >>> import icaparser as icap
        >>> extra_transcript_filters = [
                lambda x: x.get('source', '') == 'Ensembl',
                lambda x: x.get('hgnc', '') == 'KRAS'
            ]
        >>> df = pd.concat([
                     icap.get_mutation_table_for_single_file(
                         x, extra_transcript_filters=extra_transcript_filters
                     )
                for x in json_files],
               ignore_index=True)
        >>> df
    """
    # Get all positions from the json file
    positions = get_positions(file)
    # For performance reasons, we apply the extra filters first. They may, for example,
    # look for only single genes.
    # Apply extra variant filters
    for vf in extra_variant_filters:
        positions = filter_positions_by_variants(positions, vf)
    # Apply extra transcript filters
    for tf in extra_transcript_filters:
        positions = filter_positions_by_transcripts(positions, tf)
    # Keep only non-common variants
    positions = filter_positions_by_variants(
        positions,
        lambda x: (
            get_gnomad_max_af(x) <= max_af
            and get_gnomad_exome_max_af(x) <= max_af
            and get_onekg_max_af(x) <= max_af
        ),
    )
    # Keep only variants with high enough cosmic sample count
    positions = filter_positions_by_variants(
        positions, lambda x: get_cosmic_max_sample_count(x) >= min_cosmic_sample_count
    )
    # Keep only protein coding variants
    positions = filter_positions_by_transcripts(
        positions, lambda x: x.get("bioType", "") in ["protein_coding", "mRNA"]
    )
    # Keep only variants with high enough VEP impact
    positions = filter_positions_by_transcripts(
        positions,
        lambda x: get_strongest_vep_consequence_priority(x)
        <= min_vep_consequence_priority,
    )
    # Keep only canonical transcripts if requested
    if only_canonical:
        positions = filter_positions_by_transcripts(
            positions, lambda x: x.get("isCanonical", False)
        )
    # Get the transcript table
    mutation_table = get_mutation_table_for_positions(positions, hide_progress=True)
    # Done
    return mutation_table


def get_mutation_table_for_files(
    json_files: list,
    max_af: float = 0.001,
    min_vep_consequence_priority: int = 6,
    min_cosmic_sample_count: int = 0,
    only_canonical: bool = False,
    extra_variant_filters: list = [],
    extra_transcript_filters: list = [],
) -> pd.DataFrame:
    """Get an annotated table of all filtered transcripts from a list of ICA JSON files.

    Load all positions from a list of ICA JSON files and filter them.
    Positions having any remaining variants and transcripts passing the filter
    are returned as an annotated table.

    Args:
        json_files: list of ICA JSON files
        max_af: maximum allele frequency for gnomAD, gnomAD Exome and 1000 Genomes.
                Only variants with maximum allele frequencies below this threshold
                will be returned.
        min_vep_consequence_priority: only transcripts with a minimum VEP
                consequence priority not larger than this threshold will be retained.
                Consequences with priorities <= 6 change the protein sequence,
                consequences with priorities > 6 do not change the protein
                sequence.
        min_cosmic_sample_count: only variants with a maximum cosmic sample count
                not lower than this threshold will be retained
        only_canonical: if true, only canonical transcripts will be retained
        extra_variant_filters: any additional filters to apply to variants.
                Filters shall return True to keep a variant.
        extra_transcript_filters: any additional filters to apply to transcripts.
                Filters shall return True to keep a transcript.

    Returns:
        table of annotated mutations and affected transcripts.

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
    """
    df = pd.concat(
        [
            _get_mutation_table_for_single_file(
                x,
                max_af=max_af,
                min_vep_consequence_priority=min_vep_consequence_priority,
                min_cosmic_sample_count=min_cosmic_sample_count,
                only_canonical=only_canonical,
                extra_variant_filters=extra_variant_filters,
                extra_transcript_filters=extra_transcript_filters,
            )
            for x in tqdm(json_files)
        ],
        ignore_index=True,
    )
    # If none of the positions passes the filters, return the empty dataframe.
    # For an empty dataframe, the rest of this function would fail because the
    # empty dataframe does not contain these columns.
    if df.empty:
        return df
    # Make chromosome a categorical veriable set sorting order
    df["chromosome"] = df["chromosome"].astype("category")
    df["chromosome"] = df["chromosome"].cat.reorder_categories(
        natsorted(df["chromosome"].cat.categories), ordered=True
    )
    # Make isCanonical an ordered categorical variable
    df["isCanonical"] = pd.Categorical(df["isCanonical"], ordered=True)
    return df


def explode_consequence(
    mutation_table: pd.DataFrame, inplace: bool = False
) -> pd.DataFrame:
    """Explode the VEP consequence column of a mutation table.

    Exploding the VEP consequence column with the standard Pandas `explode()`
    function would return consquences as strings, not as ordered categories.
    This function will instead return a consequence column which is an ordered
    category. The categories are ordered by their impact.

    Exploding means that if a row of the input table has multiple consequences
    in the consequence column, the list of consequences will be split into
    single consequences and the output table will have multiple rows with a
    single consequence per row.

    Args:
        mutation_table: the mutation table to explode
        inplace: if True, then modify the mutation_table in place instead
            of returning a new object

    Returns:
        new mutation table with exploded consequences.

    Examples:
        >>> import icaparser as icap
        >>> icap.explode_consequence(mutation_table, inplace=True)

        >>> mutation_table_exploded = icap.explode_consequence(mutation_table)
    """
    vep_csq_cats = list(vep_csq.index)
    vep_csq_cats.reverse()
    if not inplace:
        mutation_table = deepcopy(mutation_table)
    mutation_table = mutation_table.explode(column="consequence")
    mutation_table["consequence"] = pd.Categorical(
        mutation_table["consequence"], categories=vep_csq_cats, ordered=True
    )
    return mutation_table


def add_gene_types(positions: list) -> list:
    """Adds the gene type to each transcript.

    Transcripts will be annotated with the gene type (*oncogene*, *tsg*,
    *mixed*) by adding a new attribute `geneType`. Only transcripts with
    one of these three gene types get this additional annotation. Other
    transcripts will not get the `geneType` attribute.

    Args:
        positions: list of filtered or unfiltered positions from JSON files.

    Returns:
        list of positions with additional annotation of transcripts.

    Examples:
        >>> import icaparser as icap
        >>> positions = icap.add_gene_types(positions)
    """
    for position in positions:
        for variant in position.get("variants", []):
            for transcript in variant.get("transcripts", []):
                gene = transcript.get("hgnc", "")
                gene_type = get_gene_type(gene)
                if gene_type != "":
                    transcript["geneType"] = gene_type
    return positions


def get_default_gene_type_map() -> dict:
    """Returns the default gene type map.

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
        mappings from gene types to canonical gene types.

    Examples:
        >>> import icaparser as icap
        >>> icap.get_default_gene_type_map()

    """
    # Define some mappings of gene types to rule set keys
    gene_type_map = {
        "oncogene": {"gof"},
        "tsg": {"lof"},
        "mixed": {"gof", "lof"},
        "Act": {"gof"},
        "LoF": {"lof"},
        "ambiguous": {"gof", "lof"},
        "": {"gof", "lof"},
    }
    return gene_type_map


def get_default_mutation_classification_rules(cosmic_threshold: int = 10) -> dict:
    """Returns the default rules for classifying mutations.

     Defines the default rules for classifying mutations. The returned dictionary
     has keys "gof" and "lof", and the respective values are the rule sets for
     these gene types. Each rule set is a dictionary with the keys "mutated" and
     "uncertain". The values for "mutated" or "uncertain" are dictionaries with
     three filter functions, a "position_filter", a "variant_filter", and a
     "transcript_filter". For example, a transcript will be called "mutated" if
     all three filters for "mutated" return True, and it will be called
     "uncertain", if all three filter functions for "uncertain" return True.

     These are the default rules returned by this function:

    **GOF**

    *mutated:* non-deleterious hotspot mutations.

    - position_filter: all positions retained (no restrictions by position).
    - variant_filter: keep only hotspot variants with a Cosmic sample count >=
        `cosmic_threshold`.
    - transcript_filter: keep only amino acid sequence modifying variants that
        are not most likely deleterious. This includes missense mutations and
        in-frame insertions and deletions.

    *uncertain:* non-deleterious mutations that aren't hotspots.

    - position_filter: all positions retained (no restrictions by position).
    - variant_filter: keep only non-hotspot variants with a Cosmic sample count <
        `cosmic_threshold`.
    - transcript_filter: keep only amino acid sequence modifying variants that
        are not most likely deleterious. This includes missense mutations and
        in-frame insertions and deletions.

    **LOF**

    *mutated:* deleterious mutations (such as truncations, start or stop codon
    loss).

    - position_filter: all positions retained (no restrictions by position).
    - variant_filter: all variants retained (no restrictions by variant).
    - transcript_filter: keep only deleterious variants, such as
        truncations or stop codon loss or start codon loss.

    *uncertain:* amino acid sequence modifying mutations that are not most
    likely deleterious. This includes missense mutations and in-frame insertions
    and deletions.

    - position_filter: all positions retained (no restrictions by position).
    - variant_filter: all variants retained (no restrictions by variant).
    - transcript_filter: keep only amino acid sequence modifying variants that
        are not most likely deleterious. This includes missense mutations and
        in-frame insertions and deletions.

    Args:
        cosmic_threshold: for "gof" genes, this is the "hotspot threshold" for
            Cosmic, i.e., the minimum number of samples in Cosmic having that
            mutation to consider a mutation a hot spot and, therefore, call the
            mutation "mutated". If the number of Cosmic samples is smaller, the
            mutation is called "uncertain".

    Returns:
        default mutation classification rules.

    Examples:
        >>> import icaparser as icap
        >>> icap.get_default_mutation_classification_rules()
        >>> icap.get_default_mutation_classification_rules(cosmic_threshold=20)
    """
    classification_rules = {
        "gof": {
            # potentially activating mutations (inframe indels, missense mutations)
            # that are hotspots. Currently, hotspots are defined by Cosmic only.
            "mutated": {
                "position_filter": lambda p: True,
                "variant_filter": lambda v: (
                    get_cosmic_max_sample_count(v) >= cosmic_threshold
                ),
                "transcript_filter": lambda t: (
                    get_strongest_vep_consequence_priority(t) > 4
                    and get_strongest_vep_consequence_priority(t) <= 6
                ),
            },
            # potentially activating mutations (inframe indels, missense mutations)
            # that are not hotspots
            "uncertain": {
                "position_filter": lambda p: True,
                "variant_filter": lambda v: (
                    get_cosmic_max_sample_count(v) < cosmic_threshold
                ),
                "transcript_filter": lambda t: (
                    get_strongest_vep_consequence_priority(t) > 4
                    and get_strongest_vep_consequence_priority(t) <= 6
                ),
            },
        },
        "lof": {
            # disruptive mutations (transcript_ablation, splice_acceptor_variant,
            # splice_donor_variant, stop_gained, stop_lost, start_lost,
            # frameshift_variant)
            "mutated": {
                "position_filter": lambda p: True,
                "variant_filter": lambda v: True,
                "transcript_filter": lambda t: (
                    get_strongest_vep_consequence_priority(t) <= 4
                ),
            },
            # potentially disruptive mutations (inframe indels, missense mutations)
            "uncertain": {
                "position_filter": lambda p: True,
                "variant_filter": lambda v: True,
                "transcript_filter": lambda t: (
                    get_strongest_vep_consequence_priority(t) > 4
                    and get_strongest_vep_consequence_priority(t) <= 6
                ),
            },
        },
    }
    return classification_rules


def apply_mutation_classification_rules(
    positions: list,
    rule_set: dict = get_default_mutation_classification_rules(),
    gene_type_map: dict = get_default_gene_type_map(),
    hide_progress: bool = False,
) -> tuple[list, dict]:
    """Applies mutation classification rules to all positions.

    Each variant is categorized for each transcript that overlaps with the
    genomic position of the variant. Each transcript that passes the "mutated"
    or "uncertain" mutation classification rules gets a new attribute
    `mutation_status` with the value "mutated" or "uncertain". The input list of
    positions is modified by adding the `mutation_status` attribute to
    transcripts, and the modified list of positions is returned as the first
    element of the returned tuple.

    In addition to modifying and returning the list of positions, this function
    also returns the assembled mutation status after aggregating the impact on
    all transcripts covering a variant. This is returned as the second item of
    the returned tuple. The impact depends on the type of gene ("gof" or "lof"),
    so the impacts are assembled separately for each gene type.

    The impact of a particular mutational variant can be different for different
    overlapping transcript variants of a gene, and the transcript variants can
    also belong to different genes. The strongest impact on any overlapping
    transcript of a gene is defined as the impact of that mutational variant on
    the gene. The analyst must decide which isoforms are used to classify genes.
    For example, only canonical transcripts may be considered. Alternatively,
    all transcripts or a subset of transcripts may be used. Therefore, it is
    necessary to first apply transcript-level filters to all genomic positions
    before this function is called for determining the mutation status of genes.

    The returned value is a multi-dimensional dictionary:

    `sample_id` → `gene` → `gene_type` → `variant_id` → `mutationStatus`

    Args:
        positions: list of positions.

        rule_set: rules for classifying "gof" and "lof" genes. See the default
            value for an example if a custom rule set is needed.
        gene_type_map: dictionary for mapping gene types to canonical gene
            types. See the default value for an example if a custom rule set is
            needed.

    Returns:
        A list of positions and a dictionary with assembled and aggregated
            mutations.

    Examples:
        >>> import icaparser as icap
        >>> positions, sample_muts = icap.apply_mutation_classification_rules(positions)
    """
    sample_muts = {}
    mut_rank = {"mutated": 0, "uncertain": 1, "not_mutated": 2, "": 3}
    for gene_type, rules in tqdm(
        rule_set.items(), disable=hide_progress, desc="Level 1: Gene type"
    ):
        for mut_status, filters in tqdm(
            rules.items(), disable=hide_progress, desc="Level 2: Mutation status"
        ):
            for p in tqdm(positions, disable=hide_progress, desc="Level 3: Positions"):
                if not filters["position_filter"](p):
                    continue
                sample_ids = [x["sampleId"] for x in p["samples"]]
                for v in p.get("variants", []):
                    if not filters["variant_filter"](v):
                        continue
                    for t in v.get("transcripts", []):
                        if not filters["transcript_filter"](t):
                            continue
                        # Get the gene symbol
                        gene = t.get("hgnc", "")
                        # Get the gene type of the current transcript or set it to ambiguous
                        # if the transcript has no gene type annotation
                        this_gene_type = t.get("geneType", "ambiguous")
                        # Map the gene type to {'gof'}, {'lof'} or {'gof', 'lof'}
                        this_gene_type = gene_type_map.get(
                            this_gene_type, {"gof", "lof"}
                        )
                        # Skip this transcript if the gene type for the current rule is not within
                        # the gene types of the current transcript
                        if gene_type not in this_gene_type:
                            continue
                        current_mut_status = t.get("mutationStatus", "")
                        if mut_rank[mut_status] <= mut_rank[current_mut_status]:
                            t["mutationStatus"] = mut_status
                        for sample_id in sample_ids:
                            this_sample_muts = sample_muts.get(sample_id, {})
                            this_gene_muts = this_sample_muts.get(gene, {})
                            this_gene_type_mut_status = this_gene_muts.get(
                                gene_type, {}
                            )
                            this_variant_mut_status = this_gene_type_mut_status.get(
                                v["vid"], ""
                            )
                            if (
                                mut_rank[mut_status]
                                <= mut_rank[this_variant_mut_status]
                            ):
                                this_variant_mut_status = mut_status
                            this_gene_type_mut_status[v["vid"]] = (
                                this_variant_mut_status
                            )
                            this_gene_muts[gene_type] = this_gene_type_mut_status
                            this_sample_muts[gene] = this_gene_muts
                            sample_muts[sample_id] = this_sample_muts
    return positions, sample_muts


def get_default_mutation_aggregation_rules() -> dict:
    """Returns the default mutation aggregation rules.

    Two types of the mutation status of a gene are defined - allele
    level and gene level:

    - For _gof_ genes (like oncogenes) it is sufficient that one of the
      alleles of one of the relevant isoforms has an activating mutation.
    - For _lof_ genes (like tumor suppressor genes) all alleles of all
      relevant isoforms need to be functionally disrupted, either by mutations
      or by other means.
    - For _ambiguous_ or _other_ genes, the impact of a mutation is defined as
      the highest impact according to _gof_ rules and _lof_ rules.

    For _gain of function_ (_gof_) genes, the classifications at both the allele
    and gene levels are identical unless there is supplementary information
    about activating modifications beyond mutations. In contrast, for _loss of
    function_ (_lof_) genes, classifications at the allele and gene levels may
    diverge. For instance, a truncating mutation in a tumor suppressor gene
    typically disrupts the function of the affected allele. However, other
    alleles of the same gene may remain functionally active, meaning the gene as
    a whole can still be operational, unless the mutated allele is a dominant
    negative variant. For a gene to be considered completely dysfunctional, all
    its alleles must be impaired, either through additional mutations or other
    mechanisms such as copy number deletions or hypermethylation. Consequently,
    a single variant that disrupts function at the allele level does not
    necessarily imply disruption at the gene level.

    For _loss of function_ (_lof_) genes, the available information often falls
    short of allowing a reliable estimation of functional effects. As a result,
    heuristic rules must be employed, and the analyst is tasked with deciding
    whether to utilize allele-level or gene-level classifications. A _lof_ gene
    is classified as functionally disrupted at gene level (strong impact) if it
    harbors at least two mutations, each either of strong impact or of uncertain
    impact. Should a _lof_ gene possess only one such mutation, it is classified
    as having an uncertain impact at the gene level, regardless of whether the
    mutation exhibits a strong impact at the allele level. By differentiating
    the effects at both the allele and gene levels, we maintain the flexibility
    to determine in subsequent analyses how to consolidate these categories for
    further statistical evaluations.

    The function returns a dictionary containing two keys: _gof_ and
    _lof_. Associated with each key is a function that accepts a
    dictionary of counts as its input and outputs a tuple comprising two
    elements: the mutation status at the allele level and at the gene
    level. The input dictionary of counts is expected to have two keys,
    _mutated_ and _uncertain_. The value for each key represents the
    number of variants within a gene classified as _mutated_ or
    _uncertain_, respectively.

    Returns:
        the _gof_ and _lof_ allele level and gene level aggregation rules.

    Examples:
        >>> import icaparser as icap
        >>> icap.get_default_mutation_aggregation_rules()
    """

    def _aggregation_rules_gof(mut_counts):
        if mut_counts.get("mutated", 0) >= 1:
            allele_level = gene_level = "mutated"
        elif mut_counts.get("uncertain", 0) >= 1:
            allele_level = gene_level = "uncertain"
        else:
            allele_level = gene_level = "not_mutated"
        return allele_level, gene_level

    def _aggregation_rules_lof(mut_counts):
        # allele level
        if mut_counts.get("mutated", 0) >= 1:
            allele_level = "mutated"
        elif mut_counts.get("uncertain", 0) >= 1:
            allele_level = "uncertain"
        else:
            allele_level = "not_mutated"
        # gene level
        if mut_counts.get("mutated", 0) + mut_counts.get("uncertain", 0) >= 2:
            gene_level = "mutated"
        elif mut_counts.get("mutated", 0) + mut_counts.get("uncertain", 0) == 1:
            gene_level = "uncertain"
        else:
            gene_level = "not_mutated"
        return allele_level, gene_level

    mutation_aggregation_rules = {
        "gof": _aggregation_rules_gof,
        "lof": _aggregation_rules_lof,
    }

    return mutation_aggregation_rules


def get_aggregated_mutation_table(
    positions: list,
    sample_muts: dict = None,
    mutation_classification_rules: dict = get_default_mutation_classification_rules(),
    mutation_aggregation_rules: dict = get_default_mutation_aggregation_rules(),
    gene_type_map: dict = get_default_gene_type_map(),
    hide_progress: bool = False,
) -> pd.DataFrame:
    """Returns a sample-gene-mutationStatus table.

    This function applies mutation classification rules to all mutational
    variants and aggregates the mutations according to the aggregation rules.
    This results in a table with one row for each sample-gene pair. The table
    contains several columns with impacts according to _lof_ and _gof_ rules on
    allele level and gene level and with one additional column with the maximum
    impact for both allele and gene level.

    Args:
        positions: list of positions. If sample_muts is also specified, it is
            assumed that the positions have already been processed previously by
            `apply_mutation_classification_rules` and we do not have to run
            mutation classification again.
        sample_muts: if `apply_mutation_classification_rules` has been run
            before, you can use the second return value of that function as the
            sample_muts argument. This is helpful for very large datasets
            because otherwise `apply_mutation_classification_rules` will be run
            again as an internal call within `get_aggregated_mutation_table`,
            which is time consuming for very large data sets. This also means
            that if `sample_muts` is provided as an argument, the
            `mutation_classification_rules` argument is ignored and has no
            effect.
        mutation_classification_rules: rules for classifying single mutations.
            See `get_default_mutation_classification_rules()` for details.
        mutation_aggregation_rules: rules for aggregation mutations.
            See `get_default_mutation_aggregation_reles()` for details.
        gene_type_map: dictonary for mapping gene types to canonical gene types.
            See `get_default_gene_type_map()` for details.

    Returns:
        mutation table.
    """
    gene_type_cache = {}
    for p in positions:
        for v in p.get("variants", []):
            for t in v.get("transcripts", []):
                gene_type_cache[t.get("hgnc", "")] = t.get("geneType", "")

    def get_max_mut_status(mut_status):
        if "mutated" in mut_status:
            return "mutated"
        elif "uncertain" in mut_status:
            return "uncertain"
        else:
            return "not_mutated"

    if sample_muts is None:
        p, sample_muts = apply_mutation_classification_rules(
            positions,
            rule_set=mutation_classification_rules,
            gene_type_map=gene_type_map,
            hide_progress=hide_progress,
        )
    rows = []
    for sample_id, gene_muts in tqdm(
        sample_muts.items(), disable=hide_progress, desc="Samples"
    ):
        for gene, gene_type_muts in gene_muts.items():
            gene_mut_counts = {}
            for gene_type, variant_muts in gene_type_muts.items():
                gene_type_mut_counts = gene_mut_counts.get(gene_type, {})
                for vid, mut_status in variant_muts.items():
                    mut_counts = gene_type_mut_counts.get(mut_status, 0)
                    mut_counts += 1
                    gene_type_mut_counts[mut_status] = mut_counts
                gene_mut_counts[gene_type] = gene_type_mut_counts
            gene_type = gene_type_cache.get(gene)
            canonical_gene_type = gene_type_map[gene_type]
            gof_allele_status = (
                mutation_aggregation_rules["gof"](gene_mut_counts["gof"])[0]
                if "gof" in gene_mut_counts
                else ""
            )
            lof_allele_status = (
                mutation_aggregation_rules["lof"](gene_mut_counts["lof"])[0]
                if "lof" in gene_mut_counts
                else ""
            )
            gof_gene_status = (
                mutation_aggregation_rules["gof"](gene_mut_counts["gof"])[1]
                if "gof" in gene_mut_counts
                else ""
            )
            lof_gene_status = (
                mutation_aggregation_rules["lof"](gene_mut_counts["lof"])[1]
                if "lof" in gene_mut_counts
                else ""
            )
            row = {
                "sample_id": [sample_id],
                "gene": [gene],
                "geneType": [gene_type],
                "canonicalGeneType": [canonical_gene_type],
                "gof_mutated_count": [
                    gene_mut_counts.get("gof", {}).get("mutated", 0)
                    if ({"gof", ""} & canonical_gene_type)
                    else np.NaN
                ],
                "gof_uncertain_count": [
                    gene_mut_counts.get("gof", {}).get("uncertain", 0)
                    if ({"gof", ""} & canonical_gene_type)
                    else np.NaN
                ],
                "lof_mutated_count": [
                    gene_mut_counts.get("lof", {}).get("mutated", 0)
                    if ({"lof", ""} & canonical_gene_type)
                    else np.NaN
                ],
                "lof_uncertain_count": [
                    gene_mut_counts.get("lof", {}).get("uncertain", 0)
                    if ({"lof", ""} & canonical_gene_type)
                    else np.NaN
                ],
                "gof_allele_status": [gof_allele_status],
                "lof_allele_status": [lof_allele_status],
                "gof_gene_status": [gof_gene_status],
                "lof_gene_status": [lof_gene_status],
                "allele_status": [
                    get_max_mut_status((gof_allele_status, lof_allele_status))
                ],
                "gene_status": [get_max_mut_status((gof_gene_status, lof_gene_status))],
                "gene_or_allele_status": [
                    get_max_mut_status(
                        (
                            gof_allele_status,
                            lof_allele_status,
                            gof_gene_status,
                            lof_gene_status,
                        )
                    )
                ],
            }
            rows.append(row)
    agg_mutation_table = pd.DataFrame(rows)
    count_cols = [x for x in agg_mutation_table.columns if "count" in x]
    for col in count_cols:
        agg_mutation_table[col] = agg_mutation_table[col].astype(pd.Int64Dtype())
    return agg_mutation_table


def split_multi_sample_json_file(json_file: str, output_dir: str) -> None:
    """Splits a multi-sample JSON file into sample specific JSON files.

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
        json_file: the multi-sample json input file.
        output_dir: the directory where to write the single sample JSON files.
            The directory will be created if it does not exist.

    Returns:
        None.
    """
    os.makedirs(output_dir, exist_ok=True)
    is_header_line = True
    is_position_line = False
    is_first_position_line = {}
    is_gene_line = False
    is_first_gene_line = {}
    header_end = ',"positions":['
    gene_section_line = '],"genes":['
    end_line = "]}"
    ofiles = {}
    genes = {}
    with gzip.open(json_file, "rt") as fin:
        for line in fin:
            trim_line = line.strip()
            if is_header_line:
                is_header_line = False
                is_position_line = True
                header_json = trim_line[: -len(header_end)] + "}"
                header = json.loads(header_json)["header"]
                samples = header["samples"]
                for sample in samples:
                    is_first_position_line[sample] = True
                    s_header = header.copy()
                    s_header["samples"] = [sample]
                    s_header_json = json.dumps(
                        {"header": s_header}, separators=(",", ":")
                    )
                    s_header_line = s_header_json[:-1] + header_end
                    s_fname = os.path.join(output_dir, sample + ".json.gz")
                    ofiles[sample] = gzip.open(s_fname, "wt")
                    print(s_header_line, file=ofiles[sample])
                continue
            if trim_line == gene_section_line:
                is_position_line = False
                is_gene_line = True
                for sample in samples:
                    is_first_gene_line[sample] = True
                    print("", trim_line, sep="\n", file=ofiles[sample])
                continue
            elif trim_line == end_line:
                for sample in samples:
                    print("", trim_line, sep="\n", file=ofiles[sample])
                break
            else:
                if is_position_line:
                    ## remove the training ',' if there is
                    position_json = trim_line.rstrip(",")
                    position = json.loads(position_json)
                    p_samples = position["samples"]
                    for sample_idx, sample in enumerate(samples):
                        p_sample = p_samples[sample_idx].copy()
                        # Continue if this sample is not mutated here
                        if p_sample["genotype"] in ("0/0", "./."):
                            continue
                        # Get the variant indeces for this sample
                        variant_idxs = p_sample["genotype"]
                        variant_idxs = variant_idxs.split("/")
                        variant_idxs = set(variant_idxs) - set("0")
                        variant_idxs = map(lambda x: int(x) - 1, list(variant_idxs))
                        variant_idxs = sorted(variant_idxs)
                        # Get the new genotype name for this sample
                        s_gt = p_sample["genotype"]
                        s_gt = s_gt.split("/")
                        s_gt = [
                            int(x) - min(variant_idxs) if x != "0" else 0 for x in s_gt
                        ]
                        s_gt = "/".join(map(str, s_gt))
                        p_sample["genotype"] = s_gt
                        # Adjust the list of variantFrequencies
                        if len(p_sample.get("variantFrequencies", [])) > 0:
                            p_sample["variantFrequencies"] = [
                                p_sample["variantFrequencies"][i] for i in variant_idxs
                            ]
                        else:
                            p_sample["variantFrequencies"] = []
                        # Adjust the list of alleleDepths
                        if len(p_sample.get("alleleDepths", [])) > 0:
                            p_sample["alleleDepths"] = [p_sample["alleleDepths"][0]] + [
                                p_sample["alleleDepths"][i + 1] for i in variant_idxs
                            ]
                        else:
                            p_sample["alleleDepths"] = []
                        # Keep only alleles relevant for this sample
                        s_position = position.copy()
                        s_position["samples"] = [p_sample]
                        s_position["altAlleles"] = [
                            s_position["altAlleles"][i] for i in variant_idxs
                        ]
                        s_position["variants"] = [
                            s_position["variants"][i] for i in variant_idxs
                        ]
                        if is_first_position_line[sample]:
                            finish_prev_line = ""
                            genes[sample] = set()
                        else:
                            finish_prev_line = ",\n"
                        print(
                            finish_prev_line,
                            json.dumps(s_position, separators=(",", ":")),
                            sep="",
                            end="",
                            file=ofiles[sample],
                        )
                        for variant in position.get("variants", []):
                            for transcript in variant.get("transcripts", []):
                                if "hgnc" in transcript:
                                    genes[sample].add(transcript["hgnc"])
                        is_first_position_line[sample] = False
                if is_gene_line:
                    ## remove the trailing "," if there is
                    gene_json = trim_line.rstrip(",")
                    gene = json.loads(gene_json)["name"]
                    for sample_idx, sample in enumerate(samples):
                        if sample in genes and gene in genes[sample]:
                            if is_first_gene_line[sample]:
                                finish_prv_line = ""
                            else:
                                finish_prv_line = ",\n"
                            print(
                                finish_prv_line,
                                gene_json,
                                sep="",
                                end="",
                                file=ofiles[sample],
                            )
                            is_first_gene_line[sample] = False
    # Close the output files
    for ofile in ofiles.values():
        ofile.close()
