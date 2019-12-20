rule cuffquant:
    """
    Created:
        2017-09-18 15:18:04
    Aim:
        Quantifying gene and transcript expression in RNA-Seq samples can be computationally expensive. Cuffquant allows you to compute the gene and transcript expression profiles and save these profiles to files that you can analyze later with Cuffdiff or Cuffnorm. This can help you distribute your computational load over a cluster and is recommended for analyses involving more than a handful of libraries.
    Test:
        out/cuffquant/g-GRCm38_lt-fr-firststrand/star/pe_GRCm38/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1/transcripts.gtf
        out/cufflinks/g-GRCm38_lt-fr-firststrand/star/pe_GRCm38_outFilterMultimapNmax-1000/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-P-H2AL2-WT-Rep1/transcripts.gtf
    """
    input:
        bam="out/{filler}.bam",
        gtf="out/ln/annotation_ensembl/{index}.gtf"
    output:
        cxb="out/cuffquant/g-{index}_lt-{libraryType}/{filler}/abundances.cxb"
    log:
            "out/cuffquant/g-{index}_lt-{libraryType}/{filler}/log.txt"
    params:
        outdir="out/cuffquant/g-{index}_lt-{libraryType}/{filler}"
    wildcard_constraints:
        libraryType="fr-unstranded|fr-firststrand|ff-firststrand|ff-secondstrand|fr-secondstrand|ff-unstranded|transfrags"
    threads:
        MAX_THREADS
    conda:
        "../envs/cufflinks.yaml"
    shell:
        """
        cuffquant \
            --library-type {wildcards.libraryType} \
            -p {threads}  \
            -o {params.outdir} \
            {input.gtf} \
            {input.bam} 2> {log}
        """
