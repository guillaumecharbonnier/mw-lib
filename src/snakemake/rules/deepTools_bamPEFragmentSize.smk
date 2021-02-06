rule deepTools_bamPEFragmentSize:
    """
    Created:
        2016-04-08 16h04
    Modified:
        2018-06-27 12:27:16 - Updated to current filepaths.
        2019-01-22 13:27:45 - Updated with conda and extra wildcard but not tested.
    Aim:
        This tool calculates the fragment sizes for read pairs given a BAM file from paired-end sequencing.Several regions are sampled depending on the size of the genome and number of processors to estimate thesummary statistics on the fragment lengths. Properly paired reads are preferred for computation, i.e., it will only use discordant pairs if no concordant alignments overlap with a given region. The default setting simply prints the summary statistics to the screen.
    Test:
        out/deepTools/bamPEFragmentSize/Blueprint-thymic-populations-H3K27ac-merged.pdf
    """
    input:
        bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id]),
        bai = lambda wildcards: [path + '.bai' for path in eval(mwconf['ids'][wildcards.bam_list_id])]
    output:
        pdf="out/{tool}{extra}/{bam_list_id}.pdf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="deepTools/bamPEFragmentSize"
    conda:
        "../envs/deeptools.yaml"
    threads:
        16
    shell:
        "bamPEFragmentSize --bamfiles {input.bam} --histogram {output.pdf} {params.extra} --numberOfProcessors {threads}"
