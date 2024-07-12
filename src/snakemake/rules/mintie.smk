# rule mintie_extra:
#     """
#     Created:
#         2024-03-29 20:01:56
#     Test:
#         out/arriba/pe_fq_GRCh38viral_ENSEMBL104/agent/trim_-v2/ln/alias/sst/all_samples/fastq/113281_PICCL/fusions.tsv
#     """
#     input:
#         fq1 = "out/{filler}_1.fastq.gz",
#         fq2 = "out/{filler}_2.fastq.gz",
#         gtf = "out/mintie/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_gtf_id}.gtf",
#         fa = "out/mintie/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_fa_id}.fa"
#     output:
#         tsv = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/fusions.tsv"
#     params:
#         outdir = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}"
#     wildcard_constraints:
#         tool="arriba/pe_fq",
#         arriba_fa_id = "GRCh38viral",
#         arriba_gtf_id = "ENSEMBL104"
#     conda:
#         "../envs/arriba.yaml"
#     threads:
#         8
#     shell:
#         """
#         ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

#         # Using arriba_fa_id we can decide which suffix to use for blacklist, known_fusions and protein_domains:

#         if [[ {wildcards.arriba_fa_id} == *"GRCh38"* ]] || [[ {wildcards.arriba_fa_id} == *"hg38"* ]]; then
#             STEM="_hg38_GRCh38_v2.4.0"
#         elif [[ {wildcards.arriba_fa_id} == *"GRCh37"* ]] || [[ {wildcards.arriba_fa_id} == *"hg19"* ]]; then
#             STEM="_hg19_hs37d5_GRCh37_v2.4.0"
#         else
#             echo "Invalid arriba_fa_id"
#             exit 1
#         fi

#         cd {params.outdir}

#         run_arriba.sh \
#             {WDIR}/out/arriba/download_references_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id}/STAR_index_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id} \
#             {WDIR}/{input.gtf} \
#             {WDIR}/{input.fa} \
#             $ARRIBA_FILES/blacklist$STEM.tsv.gz \
#             $ARRIBA_FILES/known_fusions$STEM.tsv.gz \
#             $ARRIBA_FILES/protein_domains$STEM.gff3 \
#             {threads} \
#             {WDIR}/{input.fq1} \
#             {WDIR}/{input.fq2}
        
#         """

# # arriba \  
# # -x {input.bam \  
# # -g {input.gtf} \  
# # -a {input.fa} \  
# # -o {output.tsv} \  
# # -O {output.discarded} \  
# # -f blacklist

# rule arriba_download_references:
#     output:
#         gtf = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_gtf_id}.gtf",
#         fa = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_fa_id}.fa"
#     params:
#         outdir = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}"
#     conda:
#         "../envs/arriba.yaml"
#     shell:
#         """
#         cd {params.outdir}
#         $CONDA_PREFIX/var/lib/arriba/download_references.sh {wildcards.arriba_fa_id}+{wildcards.arriba_gtf_id}
#         """


rule mintie_build_reference_hg38:
    """
    Currently, only building the default hg38 reference. No need for other genomes that require different steps, see:
    """
    output:
        done = "out/mintie/build_reference_hg38.done"
    conda:
        "../envs/mintie.yaml"
    shell:
        """
        touch {output.done}
        """