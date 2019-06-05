import os

include: "modules/preprocessor/helper.snake"
PARAMS = PreprocessingHelper(config, "Pugwash Pharmacogenomics")

if len(config["LOCI"]) > 1:
    raise WorkflowError("Multiple loci not yet supported")

onsuccess: PARAMS.onsuccess()
onerror: PARAMS.onerror(log)

preprocessing_params = PARAMS

# additional helpers
include: "modules/phasing/helper.snake"
include: "modules/variant_calling/helper.snake"
include: "modules/haplotyping/helper.snake"
include: "modules/structural_variation/helper.snake"
include: "modules/variant_effects/helper.snake"

# don't submit the following rules to the cluster
localrules:
    all,
    fastq_to_fasta,
    haplotypes,
    link_sources_ph,
    link_sources_vc,
    link_sources_ht,
    link_sources_sv,
    link_sources_vep,
    matches,
    pick,
    tabulate,
    laa_whitelist,
    phasing_summary


def output_files():
    file_list = [
        preprocessing_params.outputs,
        PhasingHelper(config).outputs,
        VariantCallingHelper(config).outputs,
        HaplotypingHelper(config).outputs,
        SVHelper(config).outputs,
        VEPHelper(config).outputs,
        "summary/{}/haplotype_report.html".format(list(PARAMS.genes)[0]),
        "summary/{}/deletion_report.html".format(list(PARAMS.genes)[0]),
        "summary/{}/config_report.html".format(list(PARAMS.genes)[0]),
        expand("summary/{gene}/structure/{barcode}.html", gene=list(PARAMS.genes)[0], barcode=PARAMS.barcode_ids),
        expand("summary/{gene}/phasing/{barcode}.html", gene=list(PARAMS.genes)[0], barcode=PARAMS.barcode_ids),
        expand("summary/{gene}/phasing/{barcode}.chimera_scan.tsv", gene=list(PARAMS.genes)[0], barcode=PARAMS.barcode_ids),
    ]

    #if config.get("STAGE_PARAMS", {}).get("CCS_CHECK", False):
    #    file_list.append("summary/{}/missed_variants_known.txt".format(list(PARAMS.genes)[0])),
    #    file_list.append("summary/{}/missed_variants_novel.txt".format(list(PARAMS.genes)[0]))
    
    return file_list


# workflow outputs
rule all:
    input:
        output_files()


# -------------- rules for preprocessing workflow ---------------------
include: "modules/preprocessor/rules/source_data.snake"
include: "modules/preprocessor/rules/barcoding.snake"
include: "modules/preprocessor/rules/merge_subreadset.snake"
include: "modules/preprocessor/rules/demultiplex.snake"
include: "modules/preprocessor/rules/consolidate_xml.snake"
include: "modules/preprocessor/rules/ccs.snake"


# --------------- rules for phasing workflow --------------------------
include: "modules/phasing/rules/ccs_amplicon.snake"
include: "modules/phasing/rules/haplotypes.snake"

rule link_sources_ph:
    # link the preprocessor output to the phasing input
    input:
        subreads = "preprocessor/consolidated/{barcode}.bam",
        ccs = "preprocessor/CCS/{barcode}.bam"
    output:
        subreads = "phasing/inputs/subreads/{barcode}.bam",
        ccs = "phasing/inputs/ccs/{barcode}.bam"
    shell:
        "ln -s -r {input.subreads} {output.subreads} && "
        "ln -s -r {input.ccs} {output.ccs} "


# ------------------ rules for variant calling ------------------------
include: "modules/variant_calling/rules/gene_reference.snake"
include: "modules/variant_calling/rules/call_variants.snake"

rule link_sources_vc:
    # link the phasing output to the variant calling input
    input:
        "phasing/haplotypes/{barcode}.fasta"
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
        "phasing/haplotypes/{barcode}.fastq"
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
        phasing = expand("phasing/CCS_Amplicon/{barcodes}/phasing/{barcodes}.phasing.tsv", barcodes=PARAMS.barcode_ids),
        template = srcdir("templates/haplotypes.html"),
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


rule phasing_summary:
    input:
        "phasing/CCS_Amplicon/{barcode}/phasing/{barcode}.html"
    output:
        "summary/{gene}/phasing/{barcode}.html"
    shell:
        "cp {input} {output}"


rule chimera_summary:
    input:
        "phasing/haplotypes/{barcode}.chimera_scan.tsv"
    output:
        "summary/{gene}/phasing/{barcode}.chimera_scan.tsv"
    shell:
        "cp {input} {output}"
