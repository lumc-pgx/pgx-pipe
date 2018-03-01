"""
Take a haplotype results file plus a list of variants found by ccs_check
and identify any snps which were missed by LAA
"""
import os
import json
import locus_processing
import pandas as pd


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
    pos = int(row["Pos"]) + 1
    ref = row["RefBP"]
    alt = row["AltBP"]
    length = int(row["Length"])
    
    if row["Type"] == "SNP":
        return "{pos}{ref}>{alt}".format(pos=pos, ref=ref, alt=alt)
    elif row["Type"] == "Insertion":
        return "{start}_{end}ins{alt}".format(start=pos + 1, end=pos + 1 + 1, alt=alt)
    elif row["Type"] == "Deletion":
        if length == 1:
            return "{pos}del{ref}".format(pos=pos + 1, ref=ref)
        else:
            return "{start}_{end}del{ref}".format(start=pos + 1, end=pos+length, ref=ref)


def snp_haplotypes(snp):
    snp_id = next(s.id for s in gene.snps if s.g_notation == snp)
    haplotypes = [h.type for h in gene.haplotypes if snp_id in h.snps]
    return ";".join(haplotypes)


gene = locus_processing.load_locus_yaml(snakemake.input.gene)
known_snps = set(s.g_notation for s in gene.snps)

with open(snakemake.output[0], "w") as outfile:
    for barcode in snakemake.params.barcodes:
        alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
        ccs_var = load_ccs_csv(next(f for f in snakemake.input.ccs_check if barcode in f))

        if ccs_var is not None:
            laa_variants = set(v["g_notation"] for allele in alleles for v in allele["variants"])
            ccs_var["g_notation"] = ccs_var.apply(hgvs_snp_from_row, axis=1)
            ccs_snps = ccs_var[ccs_var["Freq"] >= snakemake.params.threshold]
            missed = (set(ccs_snps["g_notation"]) - laa_variants) & known_snps
        else:
            missed = set()

        print(barcode, file=outfile, end="\t")

        if len(missed) > 0:
            missed_variants = ccs_snps[ccs_snps["g_notation"].isin(missed)]
            annotated = []
            for index, row in missed_variants.iterrows():
                annotated.append("{0} ({1:.2f};{2})".format(row["g_notation"], row["Freq"], snp_haplotypes(row["g_notation"])))
            print(",".join(annotated), file=outfile, end="")

        print("", file=outfile)

