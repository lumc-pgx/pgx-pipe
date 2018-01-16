import os
import json
import locus_processing
from jinja2 import Environment, FileSystemLoader
import itertools
from collections import Counter
import datetime
from collections import defaultdict
from interval import interval
import csv


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
    

def load_laa_summary(filename):
    # load the chimera score for each allele from the laa amplicon summary file
    with open(filename, "r") as infile:
        laa_summary = csv.DictReader(infile)
        return {row["FastaName"]: row["ChimeraScore"] for row in laa_summary}


def load_subread_summary(filename):
    # load subread summary and determine number of molecules per allele
    with open(filename, "r") as infile:
        subread_summary = csv.DictReader(infile)
        mols = defaultdict(set)
        for row in subread_summary:
            allele = next(f for f in subread_summary.fieldnames if f != "SubreadId" and float(row[f]) > 0.5)
            mols[allele].add("/".join(row["SubreadId"].split("/")[:2]))
    return {allele_id: len(mols[allele_id]) for allele_id in mols}


gene = locus_processing.load_locus_yaml(snakemake.input.gene)


def summarize_alleles(barcode):
    alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
    vep = load_vep(next(f for f in snakemake.input.vep if barcode in f))
    last = load_last(next(f for f in snakemake.input.last if barcode in f))
    chimeras = load_laa_summary(next(f for f in snakemake.input.laa_allele_summary if barcode in f))
    molecules = load_subread_summary(next(f for f in snakemake.input.laa_subread_summary if barcode in f))
    
    for allele in alleles:
        info = {}
        info["id"] = allele["sequence_id"]
        info["assignment"] = " ".join(allele["haplotype"])
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
        info["chimera_score"] = chimeras[allele_id]
        info["molecules"] = molecules[allele_id]
        
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
