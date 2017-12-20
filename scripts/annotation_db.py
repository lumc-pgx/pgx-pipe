import gffutils

gffutils.create_db(
    snakemake.input[0],
    snakemake.output[0],
    **snakemake.params
)
