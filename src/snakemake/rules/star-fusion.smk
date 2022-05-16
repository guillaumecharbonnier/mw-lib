rule star_fusion_pe_extra:
    """
    Created:
        2022-02-10 20:01:56
    Test:
    snakemake -prk --rerun-incomplete --verbose --use-conda -j4 --use-singularity --singularity-args "-B `pwd`:/data/" out/star-fusion/pe/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60/star-fusion.fusion_predictions.tsv
    """
    input:
        fwd = "out/{filler}_1.fastq.gz",
        rev = "out/{filler}_2.fastq.gz",
        ctat_lib_done = "out/tar/xvzf/wget/https/data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/done"
        #index = lambda wildcards: mwconf['ids'][wildcards.index_id],
    output:
        bam    = "out/{tool}{extra}/{filler}/star-fusion.fusion_predictions.tsv"
    params:
        ctat_lib = "out/tar/xvzf/wget/https/data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
        outdir = "out/{tool}{extra}/{filler}",
        extra = params_extra,
        # extin = lambda wildcards: "--readFilesCommand zcat" if wildcards.extin == 'fastq.gz' else "",
        # extout = lambda wildcards: "--outSAMtype BAM SortedByCoordinate" if wildcards.extout == 'bam' else "",
        # genomedir = lambda wildcards: os.path.dirname(mwconf['ids'][wildcards.index_id])
    wildcard_constraints:
        tool="star-fusion/pe",
    # star-fusion conda env creation takes eternity... (shut down after 1 day of "Solving environment")
    # conda:
    #    "../envs/star-fusion.yaml"
    container:
        "docker://trinityctat/starfusion"
    threads:
        4
    shell:
        """
        STAR-Fusion \
        --left_fq /data/{input.fwd} \
        --right_fq /data/{input.rev} \
        --genome_lib_dir /data/{params.ctat_lib} \
        -O /data/{params.outdir} \
        --runThreadN {threads} \
        --FusionInspector validate \
        --examine_coding_effect \
        --denovo_reconstruct
        """
