# Pharmacogenomics analysis of long amplicon sequences

This workflow combines the individual modules of the [PharmacogenomicsPipe](https://git.lumc.nl/PharmacogenomicsPipe) into a single workflow.

## Overview
```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    rankdir=LR;
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "last_split", color = "0.29 0.6 0.85", style="rounded"];
	1[label = "pick", color = "0.31 0.6 0.85", style="rounded"];
	2[label = "source_data", color = "0.12 0.6 0.85", style="rounded"];
	3[label = "ccs", color = "0.33 0.6 0.85", style="rounded"];
	4[label = "tabulate", color = "0.36 0.6 0.85", style="rounded"];
	5[label = "ccs_check_summary", color = "0.00 0.6 0.85", style="rounded"];
	6[label = "merge", color = "0.57 0.6 0.85", style="rounded"];
	7[label = "consolidate", color = "0.02 0.6 0.85", style="rounded"];
	8[label = "laa", color = "0.40 0.6 0.85", style="rounded"];
	9[label = "ccs_check", color = "0.05 0.6 0.85", style="rounded"];
	10[label = "gene_reference", color = "0.07 0.6 0.85", style="rounded"];
	11[label = "vep_process", color = "0.43 0.6 0.85", style="rounded"];
	12[label = "link_sources_vep", color = "0.50 0.6 0.85", style="rounded"];
	13[label = "last_tab", color = "0.45 0.6 0.85", style="rounded"];
	14[label = "lastdb", color = "0.14 0.6 0.85", style="rounded"];
	15[label = "link_sources_sv", color = "0.10 0.6 0.85", style="rounded"];
	16[label = "lastal", color = "0.48 0.6 0.85", style="rounded"];
	17[label = "all", color = "0.17 0.6 0.85", style="rounded"];
	18[label = "demultiplex", color = "0.19 0.6 0.85", style="rounded"];
	19[label = "matches", color = "0.52 0.6 0.85", style="rounded"];
	20[label = "variants", color = "0.55 0.6 0.85", style="rounded"];
	21[label = "laa_summary", color = "0.21 0.6 0.85", style="rounded"];
	22[label = "link_sources_ht", color = "0.24 0.6 0.85", style="rounded"];
	23[label = "link_sources_vc", color = "0.26 0.6 0.85", style="rounded"];
	24[label = "barcoding", color = "0.60 0.6 0.85", style="rounded"];
	25[label = "last_region", color = "0.38 0.6 0.85", style="rounded"];
	26[label = "fastq_to_fasta", color = "0.64 0.6 0.85", style="rounded"];
	16 -> 0
	19 -> 1
	7 -> 3
	19 -> 4
	9 -> 5
	8 -> 5
	24 -> 6
	18 -> 7
	7 -> 8
	3 -> 9
	12 -> 11
	20 -> 12
	0 -> 13
	8 -> 15
	14 -> 16
	15 -> 16
	1 -> 17
	3 -> 17
	19 -> 17
	4 -> 17
	5 -> 17
	25 -> 17
	20 -> 17
	21 -> 17
	11 -> 17
	26 -> 17
	6 -> 18
	22 -> 19
	10 -> 20
	23 -> 20
	8 -> 21
	20 -> 22
	26 -> 23
	2 -> 24
	13 -> 25
	8 -> 26
}
```            
