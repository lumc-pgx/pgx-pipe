import os

include: "modules/preprocessor/helper.snake"
PARAMS = PreprocessingHelper(config, "Pugwash Pharmacogenomics")

if len(config["LOCI"]) > 1:
    raise WorkflowError("Multiple loci not yet supported")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror(log)

preprocessing_params = PARAMS

# additional helpers
include: "modules/variant_calling/helper.snake"
include: "modules/haplotyping/helper.snake"
include: "modules/structural_variation/helper.snake"
include: "modules/variant_effects/helper.snake"

# don't submit the following rules to the cluster
localrules:
    all,
    fastq_to_fasta,
    barcoding_summary,
    link_sources_vc,
    link_sources_ht,
    link_sources_sv,
    link_sources_vep


# workflow outputs
rule all:
    input:
        preprocessing_params.outputs,
        VariantCallingHelper(config).outputs,
        HaplotypingHelper(config).outputs,
        SVHelper(config).outputs,
        VEPHelper(config).outputs,
        "summary/{}/haplotype_report.html".format(list(PARAMS.genes)[0]),
        "summary/{}/deletion_report.html".format(list(PARAMS.genes)[0]),
        "summary/{}/config_report.html".format(list(PARAMS.genes)[0])


# -------------- rules for preprocessing workflow ---------------------
include: "modules/preprocessor/rules/source_data.snake"
include: "modules/preprocessor/rules/barcoding.snake"
include: "modules/preprocessor/rules/barcoding_summary.snake"
include: "modules/preprocessor/rules/merge_subreadset.snake"
include: "modules/preprocessor/rules/demultiplex.snake"
include: "modules/preprocessor/rules/consolidate_xml.snake"
include: "modules/preprocessor/rules/laa.snake"
include: "modules/preprocessor/rules/laa_summary.snake"
include: "modules/preprocessor/rules/laa_whitelist.snake"
include: "modules/preprocessor/rules/ccs.snake"
include: "modules/preprocessor/rules/ccs_check.snake"
include: "modules/preprocessor/rules/fastq_to_fasta.snake"
include: "modules/preprocessor/rules/ccs_check_summary.snake"


# ------------------ rules for variant calling ------------------------
include: "modules/variant_calling/rules/gene_reference.snake"
include: "modules/variant_calling/rules/call_variants.snake"

rule link_sources_vc:
    # link the LAA output to the variant calling input
    input:
        "preprocessor/LAA/{barcode}.fasta"
    output:
        "variant_calling/inputs/{gene}/{barcode}.fasta"
    shell:
        "ln -s -r {input} {output}"


# --------------------- rules for haplotyping -------------------------
include: "modules/haplotyping/rules/matches.snake"
include: "modules/haplotyping/rules/tabulate.snake"
include: "modules/haplotyping/rules/pick.snake"
    
rule link_sources_ht:
    # link the variant calling output to the haplotyping input
    input:
        "variant_calling/variants/{gene}/{barcode}.json"
    output:
        "haplotyping/inputs/{gene}/{barcode}.json"
    shell:
        "ln -s -r {input} {output}"


# ----------------- rules for structural variation --------------------
include: "modules/structural_variation/rules/lastdb.snake"
include: "modules/structural_variation/rules/lastal.snake"
include: "modules/structural_variation/rules/last_split.snake"
include: "modules/structural_variation/rules/last_tab.snake"
include: "modules/structural_variation/rules/last_region.snake"

rule link_sources_sv:
    # link the amplicon fastq to the sv input
    input:
        "preprocessor/LAA/{barcode}.fastq"
    output:
        "structural_variation/inputs/{barcode}.fastq"
    shell:
        "ln -s -r {input} {output}"


# -------------------- rules for variant effects ----------------------
include: "modules/variant_effects/rules/do_vep.snake"

rule link_sources_vep:
    input:
        "variant_calling/variants/{gene}/{barcode}.json"
    output:
        "variant_effect/inputs/{gene}/{barcode}.json"
    shell:
        "ln -s -r {input} {output}"


# --------------------------- reporting -------------------------------
rule haplotype_summary:
    input:
        haplotypes = expand("haplotyping/haplotypes/{{gene}}/{barcodes}.haplotype.txt", barcodes=PARAMS.barcode_ids),
        matches = expand("haplotyping/matches/{{gene}}/{barcodes}.matches.json", barcodes=PARAMS.barcode_ids),
        vep = expand("variant_effect/vep/{{gene}}/{barcodes}.json", barcodes=PARAMS.barcode_ids),
        last = expand("structural_variation/last_region/{barcodes}.txt", barcodes=PARAMS.barcode_ids),
        gene = config["LOCI"][0],
        template = srcdir("templates/haplotypes.html")
    output:
        "summary/{gene}/haplotype_report.html"
    params:
        barcodes = PARAMS.barcode_ids
    conda:
        "envs/report.yaml"
    script:
        "scripts/haplotype_summary.py"

rule deletion_summary:
    input:
        last = expand("structural_variation/last_region/{barcodes}.txt", barcodes=PARAMS.barcode_ids),
        gene = config["LOCI"][0],
        template = srcdir("templates/deletions.html")
    output:
        "summary/{gene}/deletion_report.html"
    params:
        barcodes = PARAMS.barcode_ids
    conda:
        "envs/report.yaml"
    script:
        "scripts/deletion_summary.py"

rule config_summary:
    input:
        template = srcdir("templates/config.html")
    output:
        "summary/{gene}/config_report.html"
    conda:
        "envs/report.yaml"
    script:
        "scripts/config_summary.py"
