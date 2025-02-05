# bcftools liftover is actually quite new and not easily deployable. I'll wait for now...
# rule bcftools_liftover:
#     """
#     Aim:
#         Testing bcftools liftover to apply gnomAD v4 (GRCh38) as a reference for PEDIAC project
#     Test:
#         out/bcftools/liftover_chain-hg19-to-hg38/bcftools/setGT_somatic_to_missing/bcftools/fill-tags_-t_VAF/ln/updir/mw-tall-pediac-data/inp/vcf/3186patients_complet.sorted.vcf.gz
#     """
#     input:
#         vcf  = "out/{filler}.vcf.gz",
#         chain = lambda wildcards: eval(mwconf['ids'][wildcards.chain_id])
#     output:
#         vcf = "out/bcftools/liftover_{chain_id}/{filler}.vcf.gz",
#     log:
#         "out/bcftools/liftover_{chain_id}/{filler}.log"
#     benchmark:
#         "out/bcftools/liftover_{chain_id}/{filler}.benchmark.tsv"
#     conda:
#         "../../../../mw-lib/src/snakemake/envs/bcftools.yaml"
#     shell:
#         """
#         bcftools liftover --rename-chrs {input.chr_map} {input.vcf} -Oz -o {output.vcf} &> {log}
#         tabix -p vcf {output.vcf}
#         """

