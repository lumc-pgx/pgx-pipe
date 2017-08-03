# load config file
configfile: srcdir("config.yaml")

# imports
import os
import yaml
import datetime

# yaml representer for dumping config
from yaml.representer import Representer
import collections


# handlers for workflow exit status
onsuccess:
    print("Pugwash workflow completed successfully")
    yaml.add_representer(collections.OrderedDict, Representer.represent_dict)
    config_file = "config.{}.yaml".format("{:%Y-%m-%d_%H:%M:%S}".format(datetime.datetime.now()))
    with open(config_file, "w") as outfile:
        print(yaml.dump(config, default_flow_style=False), file=outfile)

onerror:
    print("Error encountered while executing workflow")
    shell("cat {log}")


include: "modules/preprocessor/helper.snake"
preproc = Preprocessing(config, "preprocessing") 

include: "modules/variant_calling/helper.snake"
var_call = VariantCalling(config, "variant calling")

include: "modules/haplotyping/helper.snake"
haplotyper = Haplotyping(config, "haplotyping")


# main workflow
localrules:
    all, preprocessing, variant_calling, haplotyping


rule all:
    input:
        expand("haplotyping/{targets}", targets=haplotyper.outputs)


rule preprocessing:
    input:
        preproc.inputs
    output:
        expand("preprocessor/{targets}", targets=preproc.outputs)
    params:
        mod_dir = srcdir("modules/preprocessor"),
        config = srcdir("config.yaml"),
        out_dir = os.path.join(os.getcwd(), "preprocessor")
    shell:
        "cd {params.mod_dir} && "
        "pipe-runner --directory {params.out_dir} --configfile {params.config}"


rule variant_calling:
    input:
        rules.preprocessing.output
    output:
        expand("variant_calling/{targets}", targets=var_call.outputs)
    params:
        mod_dir = srcdir("modules/variant_calling"),
        config = srcdir("config.yaml"),
        out_dir = os.path.join(os.getcwd(), "variant_calling"),
        laa_dir = os.path.join(os.getcwd(), "preprocessor/LAA")
    shell:
        "cd {params.mod_dir} && "
        "pipe-runner --directory {params.out_dir} --configfile {params.config} --extraconfig ALLELE_FASTA_FOLDER={params.laa_dir}"


rule haplotyping:
    input:
        rules.variant_calling.output
    output:
        expand("haplotyping/{targets}", targets=haplotyper.outputs)
    params:
        mod_dir = srcdir("modules/haplotyping"),
        config = srcdir("config.yaml"),
        out_dir = os.path.join(os.getcwd(), "haplotyping"),
        var_dir = os.path.join(os.getcwd(), "variant_calling/variants")
    shell:
        "cd {params.mod_dir} && "
        "pipe-runner --directory {params.out_dir} --configfile {params.config} --extraconfig VARIANT_DATA_PATH={params.var_dir}"

