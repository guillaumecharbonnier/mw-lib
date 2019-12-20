rule cufflinks:
    """
    Created:
        2017-02-16 16:10:02
    Aim:
        Identify new transcripts
    Doc:
        Supported library types:
            ff-firststrand
            ff-secondstrand
            ff-unstranded
            fr-firststrand
            fr-secondstrand
            fr-unstranded (default)
            transfrags
    Test:
        out/cufflinks/g-GRCm38_lt-fr-firststrand/star/pe_GRCm38/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1/transcripts.gtf
        Warning: Why does some resulting CUFF.XXX have 0 FPKM whereas they have been de novo called by Cufflinks ??????

        out/cufflinks/g-GRCm38_lt-fr-firststrand/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1/transcripts.gtf
    """
    input:
        bam="out/{filler}.bam",
        gtfGuide="out/ln/annotation_ensembl/{index}.gtf"
    output:
        gtf="out/cufflinks/g-{index}_lt-{libraryType}/{filler}/transcripts.gtf"
    log:
            "out/cufflinks/g-{index}_lt-{libraryType}/{filler}/log.txt"
    params:
        outdir="out/cufflinks/g-{index}_lt-{libraryType}/{filler}"
    wildcard_constraints:
        libraryType="fr-unstranded|fr-firststrand|ff-firststrand|ff-secondstrand|fr-secondstrand|ff-unstranded|transfrags"
    threads:
        MAX_THREADS
    conda:
        "../envs/cufflinks.yaml"
    shell:
        """
        cufflinks \
            --library-type {wildcards.libraryType} \
            --GTF-guide {input.gtfGuide} \
            -p {threads}  \
            -o {params.outdir} \
            {input.bam} 2> {log}
        """
