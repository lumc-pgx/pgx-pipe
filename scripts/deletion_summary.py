import os
import json
import locus_processing
from jinja2 import Environment, FileSystemLoader
import itertools
from collections import Counter
import datetime
from collections import defaultdict
from interval import interval

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
    last = load_last(next(f for f in snakemake.input.last if barcode in f))
    
    for allele in last:
        info = {}
        info["id"] = allele
        
        num_splits = len(last[allele])
        artifact = False
        disjoint = False
        gaps = []
        if num_splits > 1:
            # consider the allele an artifact if part of it aligns in one direction
            # and part aligns to the other direction
            artifact = len(set(x[8] for x in last[allele])) > 1
            intervals = interval()
            for split in last[allele]:
                intervals |= interval([int(split[6]), int(split[7])])
            disjoint = len(intervals) > 1
            
            gaps = []
            if disjoint:
                sorted_intervals = sorted(intervals)
                for i in range(1, len(sorted_intervals)):
                    end = sorted_intervals[i][0]
                    start = sorted_intervals[i-1][1]
                    diff = end - start
                    chrom = last[allele][i][4]
                    gaps.append("{}:{}-{} ({})".format(
                                 chrom, start, end, diff))
            
        info["artifact"] = "1" if artifact else "0"
        info["disjoint"] = "1" if disjoint else "0"
        info["gaps"] = ",".join(gaps)
    
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
            title="Pharmacogenomics Pipeline Deletion Report",
            gene=gene.name,
            sources = snakemake.config["SOURCE_DATA_PATHS"],
            summary=allele_summary,
            timestamp="{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now())
        ),
        file=outfile
    )
