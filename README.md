# Pharmacogenomics analysis of long amplicon sequences
This Snakemake workflow assigns fully phased PGx haplotypes to targetted amplicon sequencing data generated using the 
PacBio RS2 or Sequel platforms.


## Requirements
- [Conda/Miniconda](https://conda.io/miniconda.html)
- [git](https://git-scm.com/)


## Installation
- Clone this repository
  - `git clone https://git.lumc.nl/PharmacogenomicsPipe/pugwash.git`

- Change to the pugwash directory
  - `cd pugwash`

- Initialize submodules
  - `git submodule init`

- Prepare submodules
  - `git submodule update`

- Create a conda environment for running the pipeline
  - `conda env create -n pugwash -f environment.yaml`

- In order to use the pipeline on the cluster, update your .profile to use the drmaa library:
  - `echo "export DRMAA_LIBRARY_PATH=libdrmaa.so.1.0" >> ~/.profile`
  - `source ~/.profile`


## Configuration
The pipeline behavior is configured by editing [config.yaml](config.yaml).  


## Execution
- Activate the conda environment created during installation
  - `source activate pugwash`

- For parallel execution on the cluster, writing output to the default *output* directory
  - `pipe-runner`

- To specify that the pipeline should write output to a location other than the default:
  - `pipe-runner --directory path/to/output/directory`

- Note
  - The first run of the pipeline will create conda environments where required for the individual
    pipeline rules. This process can take some time, especially on slow filesystems.
  - Subsequent runs of the pipeline will re-use these environments, resulting in faster execution.


## Overview
Pugwash combines the individual analysis modules of the [PharmacogenomicsPipe project](https://git.lumc.nl/PharmacogenomicsPipe)
into a meta-workflow which can be run with a single command.

```plantuml
digraph snakemake_dag {
	graph [bb="0,0,1034,1116",
		bgcolor=white,
		margin=0
	];
    rankdir=LR;
	node [fontname=sans,
		fontsize=10,
		label="\N",
		penwidth=2,
		shape=box,
		style=rounded
	];
	edge [color=grey,
		penwidth=2
	];
	0	 [color="0.00 0.6 0.85",
		height=0.5,
		label=phasing_summary,
		pos="52,522",
		width=1.4514];
	12	 [color="0.44 0.6 0.85",
		height=0.5,
		label=all,
		pos="637,18",
		width=0.75];
	0 -> 12	 [pos="e,609.91,18.467 56.25,503.63 62.386,476.99 73,424.36 73,379 73,379 73,379 73,161 73,116.47 82.197,97.063 119,72 197.29,18.686 494.5,\
17.141 599.76,18.336"];
	1	 [color="0.40 0.6 0.85",
		height=0.5,
		label=vep_process,
		pos="572,162",
		width=1.0625];
	1 -> 12	 [pos="e,629.74,36.16 580.62,143.9 585.77,133.56 592.39,120.09 598,108 607.62,87.278 617.97,63.574 625.6,45.842"];
	32	 [color="0.38 0.6 0.85",
		height=0.5,
		label=haplotype_summary,
		pos="532,90",
		width=1.5903];
	1 -> 32	 [pos="e,541.77,108.1 562.11,143.7 557.51,135.64 551.94,125.89 546.85,116.98"];
	2	 [color="0.02 0.6 0.85",
		height=0.5,
		label=lastdb,
		pos="354,450",
		width=0.75];
	16	 [color="0.24 0.6 0.85",
		height=0.5,
		label=lastal,
		pos="354,378",
		width=0.75];
	2 -> 16	 [pos="e,354,396.1 354,431.7 354,423.98 354,414.71 354,406.11"];
	3	 [color="0.04 0.6 0.85",
		height=0.5,
		label=last_tab,
		pos="350,234",
		width=0.75694];
	28	 [color="0.36 0.6 0.85",
		height=0.5,
		label=last_region,
		pos="350,162",
		width=0.95139];
	3 -> 28	 [pos="e,350,180.1 350,215.7 350,207.98 350,198.71 350,190.11"];
	4	 [color="0.46 0.6 0.85",
		height=0.5,
		label=basic_annotations,
		pos="213,234",
		width=1.4444];
	7	 [color="0.32 0.6 0.85",
		height=0.5,
		label=annotation_db,
		pos="211,162",
		width=1.1944];
	4 -> 7	 [pos="e,211.49,180.1 212.51,215.7 212.29,207.98 212.02,198.71 211.77,190.11"];
	5	 [color="0.42 0.6 0.85",
		height=0.5,
		label=ccs_amplicon,
		pos="413,594",
		width=1.1181];
	5 -> 0	 [pos="e,104.58,536.71 372.4,586.61 315.52,577.44 208.97,559.5 119,540 117.53,539.68 116.05,539.36 114.56,539.02"];
	9	 [color="0.12 0.6 0.85",
		height=0.5,
		label=haplotypes,
		pos="707,522",
		width=0.95139];
	5 -> 9	 [pos="e,672.7,531.17 453.26,583.41 508.15,570.34 606.13,547.02 662.93,533.49"];
	5 -> 12	 [pos="e,609.91,20.19 372.58,590.73 295.04,583.97 133,557.15 133,451 133,451 133,451 133,161 133,120.25 124.34,98.836 155,72 188.33,42.823 \
492.71,25.795 599.88,20.663"];
	5 -> 32	 [pos="e,503.19,108.15 413,575.95 413,549.29 413,496.11 413,451 413,451 413,451 413,233 413,181.41 459.72,138.72 494.61,114.03"];
	6	 [color="0.06 0.6 0.85",
		height=0.5,
		label=variants,
		pos="597,378",
		width=0.75694];
	6 -> 12	 [pos="e,609.68,18.271 569.86,359.95 557.08,350.7 542.58,338.25 533,324 468.26,227.74 396.1,165.33 465,72 495.82,30.259 559.02,20.368 599.63,\
18.572"];
	15	 [color="0.22 0.6 0.85",
		height=0.5,
		label=link_sources_ht,
		pos="703,306",
		width=1.25];
	6 -> 15	 [pos="e,676.97,324.19 622.93,359.88 636.66,350.81 653.71,339.55 668.54,329.76"];
	17	 [color="0.51 0.6 0.85",
		height=0.5,
		label=link_sources_vep,
		pos="591,306",
		width=1.3472];
	6 -> 17	 [pos="e,592.47,324.1 595.52,359.7 594.86,351.98 594.06,342.71 593.32,334.11"];
	8	 [color="0.10 0.6 0.85",
		height=0.5,
		label=sv_visualization,
		pos="210,90",
		width=1.2569];
	7 -> 8	 [pos="e,210.24,108.1 210.75,143.7 210.64,135.98 210.51,126.71 210.39,118.11"];
	8 -> 12	 [pos="e,609.62,21.697 255.38,75.459 260.29,74.198 265.23,73.014 270,72 389.24,46.616 532.89,29.785 599.59,22.743"];
	9 -> 12	 [pos="e,664.2,23.64 729.46,503.97 758.67,479.71 806,432.19 806,379 806,379 806,379 806,161 806,117.99 798.7,102.12 768,72 742.34,46.827 \
703.15,33.108 674.31,25.978"];
	14	 [color="0.18 0.6 0.85",
		height=0.5,
		label=link_sources_sv,
		pos="487,450",
		width=1.2569];
	9 -> 14	 [pos="e,532.58,465.5 672.62,510.06 637.72,498.96 583.12,481.58 542.16,468.55"];
	18	 [color="0.26 0.6 0.85",
		height=0.5,
		label=chimera_summary,
		pos="888,378",
		width=1.4792];
	9 -> 18	 [pos="e,878.43,396.05 741.29,509.85 765.03,501 796.47,486.97 820,468 841.96,450.3 860.86,424.09 873.12,404.69"];
	31	 [color="0.65 0.6 0.85",
		height=0.5,
		label=link_sources_vc,
		pos="707,450",
		width=1.2569];
	9 -> 31	 [pos="e,707,468.1 707,503.7 707,495.98 707,486.71 707,478.11"];
	10	 [color="0.14 0.6 0.85",
		height=0.5,
		label=matches,
		pos="710,234",
		width=0.8125];
	10 -> 12	 [pos="e,664.43,23.577 722.08,215.62 742.08,184.48 777.11,118.15 748,72 731.91,46.486 699.78,32.992 674.21,26.021"];
	22	 [color="0.53 0.6 0.85",
		height=0.5,
		label=tabulate,
		pos="712,90",
		width=0.76389];
	10 -> 22	 [pos="e,711.76,108.19 710.24,215.87 710.58,191.67 711.21,147.21 711.61,118.39"];
	23	 [color="0.55 0.6 0.85",
		height=0.5,
		label=pick,
		pos="656,162",
		width=0.75];
	10 -> 23	 [pos="e,669.19,180.1 696.65,215.7 690.24,207.39 682.44,197.28 675.39,188.14"];
	11	 [color="0.16 0.6 0.85",
		height=0.5,
		label=source_data,
		pos="979,1098",
		width=1.0556];
	13	 [color="0.48 0.6 0.85",
		height=0.5,
		label=barcoding,
		pos="979,1026",
		width=0.88889];
	11 -> 13	 [pos="e,979,1044.1 979,1079.7 979,1072 979,1062.7 979,1054.1"];
	25	 [color="0.57 0.6 0.85",
		height=0.5,
		label=merge,
		pos="979,954",
		width=0.75];
	13 -> 25	 [pos="e,979,972.1 979,1007.7 979,999.98 979,990.71 979,982.11"];
	14 -> 16	 [pos="e,381.27,393.35 454.46,431.88 434.99,421.63 410.21,408.58 390.13,398.01"];
	15 -> 10	 [pos="e,708.29,252.1 704.73,287.7 705.5,279.98 706.43,270.71 707.29,262.11"];
	26	 [color="0.61 0.6 0.85",
		height=0.5,
		label=last_split,
		pos="350,306",
		width=0.8125];
	16 -> 26	 [pos="e,350.98,324.1 353.01,359.7 352.57,351.98 352.04,342.71 351.55,334.11"];
	17 -> 1	 [pos="e,574.3,180.19 588.71,287.87 585.47,263.67 579.52,219.21 575.67,190.39"];
	18 -> 12	 [pos="e,664.12,18.981 878.88,359.61 866,333.36 844,281.71 844,235 844,235 844,235 844,161 844,120.03 848.7,101.24 820,72 781.42,32.698 \
715.61,22.1 674.31,19.493"];
	19	 [color="0.28 0.6 0.85",
		height=0.5,
		label=demultiplex,
		pos="979,882",
		width=1];
	24	 [color="0.34 0.6 0.85",
		height=0.5,
		label=consolidate,
		pos="979,810",
		width=0.97917];
	19 -> 24	 [pos="e,979,828.1 979,863.7 979,855.98 979,846.71 979,838.11"];
	20	 [color="0.30 0.6 0.85",
		height=0.5,
		label=gene_reference,
		pos="597,450",
		width=1.2778];
	20 -> 6	 [pos="e,597,396.1 597,431.7 597,423.98 597,414.71 597,406.11"];
	21	 [color="0.20 0.6 0.85",
		height=0.5,
		label=link_sources_ph,
		pos="922,666",
		width=1.2778];
	21 -> 5	 [pos="e,453.25,600.54 875.93,658.66 780.56,645.55 562.18,615.52 463.29,601.92"];
	22 -> 12	 [pos="e,655.33,36.104 693.46,71.697 684.2,63.05 672.84,52.449 662.74,43.027"];
	23 -> 12	 [pos="e,639.3,36.189 653.71,143.87 650.47,119.67 644.52,75.211 640.67,46.393"];
	23 -> 32	 [pos="e,562.19,108.04 628.79,145.64 611.88,136.09 589.86,123.67 571.1,113.07"];
	24 -> 12	 [pos="e,664.06,19.459 993.16,791.85 1000.7,781.93 1009.5,768.89 1015,756 1030.9,718.8 1034,707.45 1034,667 1034,667 1034,667 1034,161 \
1034,118.81 1033.8,98.532 1001,72 951.29,31.792 756.68,22.058 674.23,19.725"];
	24 -> 21	 [pos="e,924.18,684.32 964.94,791.8 957.45,781.87 948.64,768.83 943,756 934.29,736.2 928.91,712.26 925.77,694.2"];
	27	 [color="0.59 0.6 0.85",
		height=0.5,
		label=ccs,
		pos="979,738",
		width=0.75];
	24 -> 27	 [pos="e,979,756.1 979,791.7 979,783.98 979,774.71 979,766.11"];
	25 -> 19	 [pos="e,979,900.1 979,935.7 979,927.98 979,918.71 979,910.11"];
	26 -> 3	 [pos="e,350,252.1 350,287.7 350,279.98 350,270.71 350,262.11"];
	26 -> 8	 [pos="e,230.37,108.19 336.87,287.84 329.41,277.7 320.16,264.45 313,252 286.65,206.14 293.28,187.37 263,144 255.86,133.77 246.63,123.76 \
237.88,115.25"];
	27 -> 12	 [pos="e,664.14,19.4 982.44,719.59 987.41,692.9 996,640.19 996,595 996,595 996,595 996,161 996,120.55 1006.4,99.793 977,72 933.99,31.323 \
753.82,21.863 674.62,19.665"];
	27 -> 21	 [pos="e,935.93,684.1 964.91,719.7 958.08,711.3 949.74,701.07 942.24,691.86"];
	28 -> 8	 [pos="e,244.38,108.19 315.75,143.88 296.94,134.47 273.42,122.71 253.33,112.67"];
	28 -> 12	 [pos="e,609.84,18.998 329.67,143.91 309.52,124.93 284.17,94.059 303,72 340.78,27.756 520.23,20.252 599.36,19.119"];
	30	 [color="0.08 0.6 0.85",
		height=0.5,
		label=deletion_summary,
		pos="365,90",
		width=1.4583];
	28 -> 30	 [pos="e,361.33,108.1 353.71,143.7 355.36,135.98 357.35,126.71 359.19,118.11"];
	28 -> 32	 [pos="e,487.79,108 384.66,147.67 411.09,137.5 447.98,123.32 478.23,111.68"];
	29	 [color="0.63 0.6 0.85",
		height=0.5,
		label=config_summary,
		pos="920,90",
		width=1.3333];
	29 -> 12	 [pos="e,664.29,25.246 871.98,75.668 867.26,74.408 862.54,73.167 858,72 793.61,55.456 718.25,37.731 674.24,27.546"];
	30 -> 12	 [pos="e,609.68,26.031 417.77,75.42 471.31,61.641 552.82,40.665 599.9,28.549"];
	31 -> 6	 [pos="e,624.01,396.19 680.09,431.88 665.84,422.81 648.15,411.55 632.76,401.76"];
	32 -> 12	 [pos="e,611.22,36.19 557.69,71.876 571.29,62.808 588.17,51.552 602.86,41.759"];
}
```            
