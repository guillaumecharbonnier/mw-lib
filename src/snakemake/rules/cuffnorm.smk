rule cuffnorm:
    """
    Created:
        2017-09-18 16:02:28
    Aim:
        Sometimes, all you want to do is normalize the expression levels from a set of RNA-Seq libraries so that they’re all on the same scale, facilitating downstream analyses such as clustering. Expression levels reported by Cufflinks in FPKM units are usually comparable between samples, but in certain situations, applying an extra level of normalization can remove sources of bias in the data. Cuffnorm normalizes a set of samples to be on as similar scales as possible, which can improve the results you obtain with other downstream tools.
    Test:
        out/cuffnorm/g-GRCm38_lt-fr-firststrand/star/pe_GRCm38/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1/transcripts.gtf
        out/cufflinks/g-GRCm38_lt-fr-firststrand/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1/transcripts.gtf
    """
    input:
        bam="out/{filler}.bam",
        gtf="out/ln/annotation_ensembl/{index}.gtf"
    output:
        #txt="out/cuffquant/g-{index}_lt-{libraryType}/{filler}/test.txt",
    log:
            "out/cuffnorm/g-{index}_lt-{libraryType}/{filler}/log.txt"
    params:
        outdir="out/cuffnorm/g-{index}_lt-{libraryType}/{filler}"
    wildcard_constraints:
        libraryType="fr-unstranded|fr-firststrand|ff-firststrand|ff-secondstrand|fr-secondstrand|ff-unstranded|transfrags"
    threads:
        4
    conda:
        "../envs/cufflinks.yaml"
    shell:
        """
        cuffnorm \
            --library-type {wildcards.libraryType} \
            -p {threads}  \
            -o {params.outdir} \
            {input.gtf} \
            {input.bam} 2> {log}
        """
