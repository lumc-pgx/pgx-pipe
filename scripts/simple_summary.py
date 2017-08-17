import json
import locus_processing
from jinja2 import Environment
import itertools
from collections import Counter
import datetime

HTML = """
<!DOCTYPE html>
<html>
<head>
<style>
    th, td {
        text-align: center;
        padding: 2px 10px;
    }
    tr, td {
        white-space: nowrap;
    }
    .allele {
        text-align: left;
    }
</style>
    <meta charset="utf-8"/>
    <title>{{ title }}</title>
</head>
<body>
    <h1>{{ title }}</h1>
    <h2>Run data:</h2>
    <ul>
    {% for source in sources %}
        <li>{{ source }}</li>
    {% endfor %}
    </ul>
    <h2>Gene: {{ gene }}</h2>
    <h2>Allele Summary</h2>
    <table>
        <tr>
            <th>Allele</th>
            <th>Assignment</th>
            <th>Variants</th>
            <th>Significant</th>
            <th>Known</th>
            <th>Novel</th>
            <th>Sift <br>(deleterious)</th>
            <th>Polyphen <br> (*damaging)</th>
        </tr>
    {% for item in summary %}
        <tr>
            <td class='allele'>{{ item['id'] }}</td>
            <td>{{ item['assignment'] }}</td>
            <td>{{ item['total'] }}</td>
            <td>{{ item['significant'] }}</td>
            <td>{{ item['known'] }}</td>
            <td>{{ item['unknown'] }}</td>
            <td>{{ item['sift'] }}</td>
            <td>{{ item['poly'] }}</td>
        </tr>
    {% endfor %}
    </table>
    <br><br>
    <footer>
        Creation timestamp: {{ timestamp }}
    </footer>
</body>
</html>
"""


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
        
gene = locus_processing.load_locus_yaml(snakemake.input.gene)

def summarize_alleles(barcode):
    alleles = load_alleles(next(f for f in snakemake.input.haplotypes if barcode in f))
    matches = load_matches(next(f for f in snakemake.input.matches if barcode in f))
    vep = load_vep(next(f for f in snakemake.input.vep if barcode in f))
    
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
    
        yield info


allele_summary = []
for barcode in snakemake.params.barcodes:
    for allele in summarize_alleles(barcode):
        allele_summary.append(allele)
        
with open(snakemake.output[0], "w") as outfile:
    print(
        Environment().from_string(HTML).render(
            title="Pharmacogenomics Pipeline Report",
            gene=gene.name,
            sources = snakemake.config["SOURCE_DATA_PATHS"],
            summary=allele_summary,
            timestamp="{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now())
        ),
        file=outfile
    )
