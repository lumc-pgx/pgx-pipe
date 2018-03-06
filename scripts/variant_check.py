"""
Take a haplotype results file plus a list of variants found by ccs_check
and identify any snps which were missed by LAA
"""
import os
import json
import locus_processing
import pandas as pd
from variant_tools import mutalyzer, alignment, variant_calling
import pyfaidx


def load_alleles(filename):
    # load the haplotypes json file
    with open(filename, "r") as infile:
        return json.load(infile)


def load_ccs_csv(filename):
    # load the ccs_check counts table
    try:
        return pd.read_csv(filename, header=0)
    except pd.errors.EmptyDataError:
        return None


def hgvs_snp_from_row(row):
    # construct an hgvs-like representation of the ccs_check variant
    start = row["Start"]
    end = row["End"]
    ref = row["RefBP"]
    alt = row["AltBP"]
    v_type = row["Type"]
    
    if v_type == "SNP":
        return "{pos}{ref}>{alt}".format(pos=start, ref=ref, alt=alt)
    elif v_type == "Insertion":
        return "{start}_{end}ins{alt}".format(start=start, end=end, alt=alt)
    elif v_type == "Deletion":
        if start == end:
            return "{pos}del{ref}".format(pos=start, ref=ref)
        else:
            return "{start}_{end}del{ref}".format(start=start, end=end, ref=ref)


def snp_haplotypes(snp):
    # get a list of haplotypes which are associated with a specific snp
    snp_id = next(s.id for s in gene.snps if s.g_notation == snp)
    haplotypes = [h.type for h in gene.haplotypes if snp_id in h.snps]
    return ";".join(haplotypes)


def normalize_variant(row):
    # Generate a normalized representation of the ccs_check variant.
    # This uses the variant_tools module to perform local realignment of
    # indels followed by calling of variants with an hgvs-like notation
    
    # the un-normalized variant
    v = row["g_notation"]

    # do not re-normalize substitutions
    if row["Type"] == "SNP":
        return v

    # variant start/end positions
    start, end = (row["Start"], row["End"])

    # we apply the variant to 'windowsize' nt each size of the variant position
    window = snakemake.params.windowsize
   
    # obtain the sequence with the variant applied
    mutated_seq = mutalyzer.apply_variants(
        [v], 
        gene.chromosome.accession, 
        start - window, 
        end + window
    )
    
    # get the reference sequence for the same region
    ref_seq = genome[gene.chromosome.name][start - window - 1:end + window].seq
    
    # perform alignment, right-justifed
    aln = alignment.align(ref_seq, mutated_seq)
    
    # re-call the variant
    variants = variant_calling.call_variants(aln, start - window - 1)
    return str(variants[0])


# some globals
gene = locus_processing.load_locus_yaml(snakemake.input.gene)
known_snps = set(s.g_notation for s in gene.snps)
genome = pyfaidx.Fasta(snakemake.input.genome)


# do the work
with open(snakemake.output.known, "w") as known_out, \
     open(snakemake.output.novel, "w") as novel_out:
    # one barcode at a time
    for barcode in snakemake.params.barcodes:
        # list of dicts containing the allele assignments
        alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
       
        # pandas dataframe containing the ccs_check data
        ccs_var = load_ccs_csv(next(f for f in snakemake.input.ccs_check if barcode in f))

        if ccs_var is not None:
            # a set containing the variants assigned to an allele
            laa_variants = set(v["g_notation"] for allele in alleles for v in allele["variants"])
            
            # filter out any ccs_check variants < the specified frequency threshold
            ccs_var = ccs_var[ccs_var["Freq"] >= snakemake.params.threshold]

            # add columns to the dataframe containing the variant start, end and unnormalized variant representation
            ccs_var["Start"] = ccs_var.apply(lambda x: x["Pos"] + 1 if x["Type"] in ("SNP", "Insertion") else x["Pos"] + 2, axis=1).astype(int)
            ccs_var["End"] = ccs_var.apply(lambda x: x["Start"] + 1 if x["Type"] == "Insertion" else x["Start"] + x["Length"] - 1, axis=1).astype(int)
            ccs_var["g_notation"] = ccs_var.apply(hgvs_snp_from_row, axis=1)
            
            # add a column containing the normalized variant
            ccs_var["Normalized"] = ccs_var.apply(normalize_variant, axis=1)
            
            # Make a set containing the description of variants which were found by ccs_check
            # but were not included in any of the LAA sequences
            missed = set(ccs_var["Normalized"]) - laa_variants
            
            # known
            missed_known = missed & known_snps
            
            # and novel
            missed_novel = missed - known_snps
        else:
            missed_known = set()
            missed_novel = set()

        # output known
        print(barcode, file=known_out, end="\t")

        if len(missed_known) > 0:
            # filter the dataframe to include only the missed variants
            missed_variants = ccs_var[ccs_var["Normalized"].isin(missed_known)]
            # annotate and dump
            annotated = []
            for index, row in missed_variants.iterrows():
                annotated.append("{0} ({1:.2f};{2})".format(row["Normalized"], row["Freq"], snp_haplotypes(row["Normalized"])))
            print(",".join(annotated), file=known_out, end="")

        print("", file=known_out)
        
        # output novel
        print(barcode, file=novel_out, end="\t")

        if len(missed_novel) > 0:
            # filter the dataframe to include only the missed variants
            missed_variants = ccs_var[ccs_var["Normalized"].isin(missed_novel)]
            # annotate and dump
            annotated = []
            for index, row in missed_variants.iterrows():
                annotated.append("{0} ({1:.2f})".format(row["Normalized"], row["Freq"],))
            print(",".join(annotated), file=novel_out, end="")

        print("", file=novel_out)

