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
        alleles = json.load(infile)
    return alleles


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
    # generate a normalized representation of the ccs_check variant
    window = snakemake.params.windowsize
    accession = gene.chromosome.accession
    chrom = gene.chromosome.name
    
    v = row["g_notation"]
    start, end = (row["Start"], row["End"])

    mutated_seq = mutalyzer.apply_variants([v], accession, start - window, end + window)
    ref_seq = genome[chrom][start - window - 1:end + window].seq
    aln = alignment.align(ref_seq, mutated_seq)
    variants = variant_calling.call_variants(aln, start - window - 1)
    return str(variants[0])


gene = locus_processing.load_locus_yaml(snakemake.input.gene)
known_snps = set(s.g_notation for s in gene.snps)
genome = pyfaidx.Fasta(snakemake.input.genome)

with open(snakemake.output[0], "w") as outfile:
    for barcode in snakemake.params.barcodes:
        alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
        ccs_var = load_ccs_csv(next(f for f in snakemake.input.ccs_check if barcode in f))

        if ccs_var is not None:
            laa_variants = set(v["g_notation"] for allele in alleles for v in allele["variants"])
            ccs_var = ccs_var[ccs_var["Freq"] >= snakemake.params.threshold]
            ccs_var["Start"] = ccs_var.apply(lambda x: x["Pos"] + 1 if x["Type"] in ("SNP", "Insertion") else x["Pos"] + 2, axis=1).astype(int)
            ccs_var["End"] = ccs_var.apply(lambda x: x["Start"] + 1 if x["Type"] == "Insertion" else x["Start"] + x["Length"] - 1, axis=1).astype(int)
            ccs_var["g_notation"] = ccs_var.apply(hgvs_snp_from_row, axis=1)
            ccs_var["Normalized"] = ccs_var.apply(normalize_variant, axis=1)
                
            missed = (set(ccs_var["Normalized"]) - laa_variants) & known_snps
        else:
            missed = set()

        print(barcode, file=outfile, end="\t")

        if len(missed) > 0:
            missed_variants = ccs_var[ccs_var["Normalized"].isin(missed)]
            annotated = []
            for index, row in missed_variants.iterrows():
                annotated.append("{0} ({1:.2f};{2})".format(row["Normalized"], row["Freq"], snp_haplotypes(row["Normalized"])))
            print(",".join(annotated), file=outfile, end="")

        print("", file=outfile)

