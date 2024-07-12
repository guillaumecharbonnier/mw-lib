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
        tsv    = "out/{tool}{extra}/{filler}/star-fusion.fusion_predictions.tsv"
    params:
        ctat_lib = "out/tar/xvzf/wget/https/data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir",
        outdir = "out/{tool}{extra}/{filler}",
        extra = params_extra,
        # extin = lambda wildcards: "--readFilesCommand zcat" if wildcards.extin == 'fastq.gz' else "",
        # extout = lambda wildcards: "--outSAMtype BAM SortedByCoordinate" if wildcards.extout == 'bam' else "",
        # genomedir = lambda wildcards: os.path.dirname(mwconf['ids'][wildcards.index_id])
    wildcard_constraints:
        tool="star-fusion/pe",
    # Facing an error with conda during fusionInsector validate
    # Trying singularity instead. Authors recommend using Singularity
    # conda:
    #     "../envs/star-fusion.yaml"
    container:
        "docker://trinityctat/starfusion"
        # Remember to add this to the snakemake call: --use-singularity --singularity-args "--bind /mnt/thymus/synoSalva/illumina_sequencing_data/mw/mw-tall-data:/data"
    threads:
        4
    shell:
        """
        STAR-Fusion \
        --left_fq {WDIR}/{input.fwd} \
        --right_fq {WDIR}/{input.rev} \
        --genome_lib_dir {WDIR}/{params.ctat_lib} \
        -O {WDIR}/{params.outdir} \
        --FusionInspector validate \
        --examine_coding_effect \
        --denovo_reconstruct
        """
