import os
import json
import locus_processing
from jinja2 import Environment, FileSystemLoader
import itertools
from collections import Counter
import datetime
from collections import defaultdict
from interval import interval
import pandas as pd


def load_alleles(filename):
    # load the haplotypes json file
    with open(filename, "r") as infile:
        alleles = json.load(infile)
    return alleles


def load_vep(filename):
    # load the vep json file
    with open(filename, "r") as infile:
        vep = json.load(infile)
    return vep


def load_last(filename):
    # load the regions from sv module into a dict indexed by allele id
    last = defaultdict(list)
    with open(filename, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            fields = [x.strip() for x in line.split("\t")]
            last[fields[0]].append(fields[1:])
    return last
    

def load_phasing_summary(filename):
    # load phasing summary
    try:
        return pd.read_csv(filename, sep="\t")
    except pd.errors.EmptyDataError:
        return None


def count_molecules(allele_id, phasing):
    if phasing is None:
        return 0

    counts = phasing.groupby(["cluster", "phase"])["id"].count()
    fields = allele_id.split(".")
    cluster = int(fields[1].split("cluster")[-1])
    phase = int(fields[2].split("haplotype")[-1])
    return counts[(cluster, phase)]


def count_passes(allele_id, phasing):
    if phasing is None:
        return 0

    counts = phasing.groupby(["cluster", "phase"])["np"].sum()
    fields = allele_id.split(".")
    cluster = int(fields[1].split("cluster")[-1])
    phase = int(fields[2].split("haplotype")[-1])
    return counts[(cluster, phase)]


gene = locus_processing.load_locus_yaml(snakemake.input.gene)


def summarize_alleles(barcode):
    alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
    vep = load_vep(next(f for f in snakemake.input.vep if barcode in f))
    last = load_last(next(f for f in snakemake.input.last if barcode in f))
    phasing = load_phasing_summary(next(f for f in snakemake.input.phasing if barcode in f))
    
    num_alleles = len(alleles)
    first = True
    
    allele_info = []
    for allele in alleles:
        info = {}        
        info["id"] = allele["sequence_id"]
        info["assignment"] = " ".join(allele["haplotype"])
        
        info["molecules"] = count_molecules(allele["sequence_id"], phasing)
        info["passes"] = count_passes(allele["sequence_id"], phasing)
        
        match = allele["haplotypes"]
        info["known"] = len([v for v in allele["variants"] if "known" in v["tags"]])
        info["unknown"] = len([v for v in allele["variants"] if "novel" in v["tags"]])
        info["significant"] = len([v for v in allele["variants"] if "significant" in v["tags"]])
        info["total"] = info["known"] + info["unknown"]
        
        vep_effects = itertools.chain.from_iterable([v["vep"] for v in vep if v["sequence_id"] == allele["sequence_id"]])
        vep_filtered = [v for v in vep_effects if v["input"].split(".")[-1] in set(v["g_notation"] for v in allele["variants"] if "novel" in v["tags"])]
    
        consequences = [v["transcript_consequences"][0] for v in vep_filtered]
    
        sift = [v["sift_prediction"] for v in consequences if "sift_prediction" in v and v["sift_prediction"] != "tolerated"]
        poly = [v["polyphen_prediction"] for v in consequences if "polyphen_prediction" in v and "damaging" in v["polyphen_prediction"]]
    
        info["sift"] = len(sift)
        info["poly"] = len(poly)
        
        allele_id = allele["sequence_id"]
        num_splits = len(last[allele_id])
        artifact = False
        disjoint = False
        if num_splits > 1:
            # consider the allele an artifact if part of it aligns in one direction
            # and part aligns to the other direction
            artifact = len(set(x[8] for x in last[allele_id])) > 1
            intervals = interval()
            for split in last[allele_id]:
                intervals |= interval([int(split[6]), int(split[7])])
            disjoint = len(intervals) > 1

        info["artifact"] = "1" if artifact else "0"
        info["disjoint"] = "1" if disjoint else "0"

        allele_info.append(info)
    
    total_molecules = sum(i["molecules"] for i in allele_info)
    allele_info.sort(key=lambda x: x["molecules"], reverse=True)
    for info in allele_info:
        info["fraction"] = info["molecules"] / total_molecules
        yield info


allele_summary = []
shade = True;
for barcode in snakemake.params.barcodes:
    alleles = list(summarize_alleles(barcode))
    
    for allele in alleles:
        allele["shade"] = "dark" if shade else "light"
        allele_summary.append(allele)
    
    if len(alleles) > 0:
        shade = not shade


j2_env = Environment(loader=FileSystemLoader(os.path.dirname(snakemake.input.template)),
                     trim_blocks=True)

j2_template = j2_env.get_template(os.path.basename(snakemake.input.template))

with open(snakemake.output[0], "w") as outfile:
    print(
        j2_template.render(
            title="Pharmacogenomics Pipeline Report",
            gene=gene.name,
            sources = snakemake.config["SOURCE_DATA_PATHS"],
            summary=allele_summary,
            timestamp="{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now())
        ),
        file=outfile
    )
