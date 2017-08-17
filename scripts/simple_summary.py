import os
import json
import locus_processing
from jinja2 import Environment, FileSystemLoader
import itertools
from collections import Counter
import datetime
from collections import defaultdict
from interval import interval

def load_alleles(filename):
    # load the haplotype assignments into a dictionary, using allele id as key
    with open(filename, "r") as infile:
        alleles = {}
        for line in infile:
            line = [x.strip() for x in line.split("\t")]
            if len(line) > 0:
                alleles[line[0]] = " ".join(line[1:])
    return alleles

def load_matches(filename):
    # load the matches json file
    with open(filename, "r") as infile:
        matches = json.load(infile)
    return matches

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

gene = locus_processing.load_locus_yaml(snakemake.input.gene)

def summarize_alleles(barcode):
    alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
    matches = load_matches(next(f for f in snakemake.input.matches if barcode in f))
    vep = load_vep(next(f for f in snakemake.input.vep if barcode in f))
    last = load_last(next(f for f in snakemake.input.last if barcode in f))
    
    for allele in alleles:
        info = {}
        info["id"] = allele
        info["assignment"] = alleles[allele]
        match = next(m for m in matches if m["sequence_id"] == allele)
        info["known"] = len(match["known_variants"])
        info["unknown"] = len(match["novel_variants"])
        info["significant"] = len(match["significant_variants"])
        info["total"] = info["known"] + info["unknown"]
        
        vep_effects = itertools.chain.from_iterable([v["vep"] for v in vep if v["sequence_id"] == allele])
        vep_filtered = [v for v in vep_effects if v["input"].split(".")[-1] in set(match["novel_variants"])]
    
        consequences = [v["transcript_consequences"][0] for v in vep_filtered]
    
        sift = [v["sift_prediction"] for v in consequences if "sift_prediction" in v and v["sift_prediction"] != "tolerated"]
        poly = [v["polyphen_prediction"] for v in consequences if "polyphen_prediction" in v and "damaging" in v["polyphen_prediction"]]
    
        info["sift"] = len(sift)
        info["poly"] = len(poly)
        
        num_splits = len(last[allele])
        artifact = False
        disjoint = False
        if num_splits > 1:
            # consider the allele an artifact if part of it aligns in one direction
            # and part aligns to the other direction
            artifact = len(set(x[8] for x in last[allele])) > 1
            intervals = interval()
            for split in last[allele]:
                intervals |= interval([int(split[6]), int(split[7])])
            disjoint = len(intervals) > 1
            
        info["artifact"] = "1" if artifact else "0"
        info["disjoint"] = "1" if disjoint else "0"
    
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
