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
    link_sources_vc,
    link_sources_ht,
    link_sources_sv,
    link_sources_vep,
    matches,
    pick,
    tabulate,
    laa_whitelist


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
        "summary/{}/config_report.html".format(list(PARAMS.genes)[0]),
        expand("summary/{gene}/structure/{barcode}.html", gene=list(PARAMS.genes)[0], barcode=PARAMS.barcode_ids),

if config.get("STAGE_PARAMS", {}).get("CCS_CHECK", False):
    rules.all.input.append("summary/{}/missed_variants_known.txt".format(list(PARAMS.genes)[0])),
    rules.all.input.append("summary/{}/missed_variants_novel.txt".format(list(PARAMS.genes)[0]))


# -------------- rules for preprocessing workflow ---------------------
include: "modules/preprocessor/rules/source_data.snake"
include: "modules/preprocessor/rules/barcoding.snake"
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
        haplotypes = expand("haplotyping/haplotypes/{{gene}}/{barcodes}.haplotype.json", barcodes=PARAMS.barcode_ids),
        vep = expand("variant_effect/vep/{{gene}}/{barcodes}.json", barcodes=PARAMS.barcode_ids),
        last = expand("structural_variation/last_region/{barcodes}.txt", barcodes=PARAMS.barcode_ids),
        gene = config["LOCI"][0],
        template = srcdir("templates/haplotypes.html"),
        laa_allele_summary = expand("preprocessor/LAA/{barcodes}_summary.csv", barcodes=PARAMS.barcode_ids),
        laa_subread_summary = expand("preprocessor/LAA/subreads.{barcodes}--{barcodes}.csv", barcodes=PARAMS.barcode_ids),
        laa_input_summary = expand("preprocessor/LAA/{barcodes}_input.csv", barcodes=PARAMS.barcode_ids)
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

rule missed_variants:
    input:
        haplotypes = expand("haplotyping/haplotypes/{{gene}}/{barcodes}.haplotype.json", barcodes=PARAMS.barcode_ids),
        ccs_check = expand("preprocessor/summary/ccs_check/{barcodes}.csv", barcodes=PARAMS.barcode_ids),
        gene = config["LOCI"][0],
        genome = config["GENOME"]
    output:
        known = "summary/{gene}/missed_variants_known.txt",
        novel = "summary/{gene}/missed_variants_novel.txt"
    conda:
        "envs/variant_check.yaml"
    params:
        threshold=0.2,
        windowsize=50,
        barcodes = PARAMS.barcode_ids
    script:
        "scripts/variant_check.py"

rule basic_annotations:
    input:
        config["ANNOTATIONS"]
    output:
        "annotations/{}.basic.gtf".format(os.path.splitext(os.path.basename(config["ANNOTATIONS"]))[0])
    shell:
        """
        grep $'\tgene\t\|tag "basic";' {input} > {output} 
        """

rule annotation_db:
    input:
        "annotations/{}.basic.gtf".format(os.path.splitext(os.path.basename(config["ANNOTATIONS"]))[0])
    output:
        "annotations/{}.basic.gtf.sqlite".format(os.path.splitext(os.path.basename(config["ANNOTATIONS"]))[0])
    params:
        # params to pass to gffutils createdb
        # these are the recommended parameters for a gencode annotation gtf
        # see https://daler.github.io/gffutils/examples.html#gencode-v19-gtf
        keep_order=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True
    conda:
        "envs/annotations.yaml"
    script:
        "scripts/annotation_db.py"

rule sv_visualization:
    input: 
        annotations = rules.annotation_db.output[0],
        alignments = "structural_variation/last_region/{barcode}.txt",
        splits = "structural_variation/last_split/{barcode}.maf"
    output:
        "summary/{gene}/structure/{barcode}.html"
    conda:
        "envs/sv_viz.yaml"
    script:
        "scripts/sv_report.py"
        
