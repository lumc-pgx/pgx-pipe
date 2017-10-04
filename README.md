# Pharmacogenomics analysis of long amplicon sequences

This workflow combines the individual modules of the [PharmacogenomicsPipe](https://git.lumc.nl/PharmacogenomicsPipe) into a single workflow.

## Overview
```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    rankdir=LR;
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "ccs", color = "0.37 0.6 0.85", style="rounded"];
	1[label = "tabulate", color = "0.35 0.6 0.85", style="rounded"];
	2[label = "laa_summary", color = "0.06 0.6 0.85", style="rounded"];
	3[label = "last_split", color = "0.41 0.6 0.85", style="rounded"];
	4[label = "fastq_to_fasta", color = "0.02 0.6 0.85", style="rounded"];
	5[label = "consolidate", color = "0.08 0.6 0.85", style="rounded"];
	6[label = "config_summary", color = "0.45 0.6 0.85", style="rounded"];
	7[label = "deletion_summary", color = "0.47 0.6 0.85", style="rounded"];
	8[label = "lastal", color = "0.39 0.6 0.85", style="rounded"];
	9[label = "laa_1", color = "0.33 0.6 0.85", style="rounded"];
	10[label = "lastdb", color = "0.25 0.6 0.85", style="rounded"];
	11[label = "all", color = "0.10 0.6 0.85", style="rounded"];
	12[label = "ccs_check", color = "0.12 0.6 0.85", style="rounded"];
	13[label = "laa_whitelist", color = "0.61 0.6 0.85", style="rounded"];
	14[label = "link_sources_sv", color = "0.51 0.6 0.85", style="rounded"];
	15[label = "demultiplex", color = "0.14 0.6 0.85", style="rounded"];
	16[label = "last_tab", color = "0.16 0.6 0.85", style="rounded"];
	17[label = "last_region", color = "0.55 0.6 0.85", style="rounded"];
	18[label = "barcoding", color = "0.53 0.6 0.85", style="rounded"];
	19[label = "merge", color = "0.43 0.6 0.85", style="rounded"];
	20[label = "haplotype_summary", color = "0.18 0.6 0.85", style="rounded"];
	21[label = "pick", color = "0.20 0.6 0.85", style="rounded"];
	22[label = "variants", color = "0.22 0.6 0.85", style="rounded"];
	23[label = "ccs_check_summary", color = "0.57 0.6 0.85", style="rounded"];
	24[label = "source_data", color = "0.49 0.6 0.85", style="rounded"];
	25[label = "link_laa", color = "0.24 0.6 0.85", style="rounded"];
	26[label = "vep_process", color = "0.59 0.6 0.85", style="rounded"];
	27[label = "matches", color = "0.27 0.6 0.85", style="rounded"];
	28[label = "gene_reference", color = "0.29 0.6 0.85", style="rounded"];
	29[label = "laa_2", color = "0.04 0.6 0.85", style="rounded"];
	30[label = "link_sources_vc", color = "0.31 0.6 0.85", style="rounded"];
	31[label = "link_sources_ht", color = "0.63 0.6 0.85", style="rounded"];
	32[label = "link_sources_vep", color = "0.65 0.6 0.85", style="rounded"];
	5 -> 0
	27 -> 1
	25 -> 2
	8 -> 3
	25 -> 4
	15 -> 5
	17 -> 7
	14 -> 8
	10 -> 8
	5 -> 9
	0 -> 11
	17 -> 11
	4 -> 11
	20 -> 11
	21 -> 11
	22 -> 11
	23 -> 11
	2 -> 11
	7 -> 11
	6 -> 11
	26 -> 11
	27 -> 11
	1 -> 11
	0 -> 12
	9 -> 13
	25 -> 14
	19 -> 15
	3 -> 16
	16 -> 17
	24 -> 18
	18 -> 19
	17 -> 20
	21 -> 20
	26 -> 20
	27 -> 20
	27 -> 21
	30 -> 22
	28 -> 22
	25 -> 23
	12 -> 23
	29 -> 25
	32 -> 26
	31 -> 27
	13 -> 29
	5 -> 29
	4 -> 30
	22 -> 31
	22 -> 32
}
```            
