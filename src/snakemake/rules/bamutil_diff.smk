rule bamutil_diff_abara2check:
    """
    Created:
        2020-01-23 09:53:12
    Aim:
        I want to compare bam before and after abra realignement
    Note:
        "echo {input.bam}[1]; echo {input.bam}[[1]]"
    Test:
        out/bamutil/diff/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-hg19-main-chr/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.txt
    """
    input:
        bam_before="out/{filler}.bam",
        bam_after="out/samtools/index/abra2/_fa-genome-hg19-main-chr/{filler}.bam"
        #bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id])
    output:
        "out/bamutil/diff/{filler}.txt"
    conda:
        "../envs/bamutil.yaml"
    shell:
        "bam diff --in1 {input.bam_before} --in2 {input.bam_after} > {output}"
