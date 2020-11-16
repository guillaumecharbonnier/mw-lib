rule SvABA_germline:
    """
    Created:
        2020-08-18 20:54:37
    Test:
        out/SvABA_bwa-index-GRCh38/bwa/mem2_se_bwa-index-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.contigs.bam
    """
    input:
        bam="out/{filler}.bam",
        fasta = lambda wildcards: eval(config['ids'][wildcards.bwa_index_id])
    output:
        bam="out/{tool}_{bwa_index_id}/{filler}.contigs.bam",
    params:
        outdir= "out/{tool}_{bwa_index_id}/{filler}"
    wildcard_constraints:
        tool = "SvABA"
    conda:
        "../envs/svaba.yaml"
    shell:
        """
        ## Set -I to not do mate-region lookup if mates are mapped to different chromosome.
        ##   This is appropriate for germline-analysis, where we don't have a built-in control
        ##   to against mapping artifacts, and we don't want to get bogged down with mate-pair
        ##   lookups.
        ## Set -L to 6 which means that 6 or more mate reads must be clustered to 
        ##   trigger a mate lookup. This also reduces spurious lookups as above, and is more 
        ##   appropriate the expected ALT counts found in a germline sample 
        ##   (as opposed to impure, subclonal events in cancer that may have few discordant reads).
        svaba run -t {input.bam} -p {threads} -L 6 -I -a {params.outdir} -G `echo {input.fasta} | sed 's/ .*$//'`
        """

