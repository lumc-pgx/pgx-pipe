# Pharmacogenomics analysis of long amplicon sequences

This workflow combines the individual modules of the [PharmacogenomicsPipe](https://git.lumc.nl/PharmacogenomicsPipe) into a single workflow.

## Overview
```plantuml
digraph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	0[label = "variant_calling", color = "0.00 0.6 0.85", style="rounded"];
	1[label = "all", color = "0.13 0.6 0.85", style="rounded"];
	2[label = "preprocessing", color = "0.53 0.6 0.85", style="rounded"];
	3[label = "haplotyping", color = "0.40 0.6 0.85", style="rounded"];
	4[label = "structural_variation", color = "0.27 0.6 0.85", style="rounded"];
	2 -> 0
	3 -> 1
	4 -> 1
	0 -> 3
	2 -> 4
}
```            
