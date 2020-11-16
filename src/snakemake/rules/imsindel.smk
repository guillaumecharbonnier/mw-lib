rule imsindel:
    """
    Test:
        out/imsindel/_fa-genome-GRCh38-r94-chr1/samtools/index/samtools/sort/samtools/view_sam_to_bam/bwa/mem_se_fa-genome-GRCh38-r94-chr1/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.vcf
    """
    input:
        #imsindel="../IMSindel/bin/imsindel",
        bam="out/{filler}.bam",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}/1.out"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="imsindel/"
    threads:
        MAX_THREADS
    conda:
        "../envs/imsindel.yaml"
    shell:
        """
        imsindel --thread {threads} --bam {input.bam} --chr 1 --indelsize 10000 --outd `dirname {output}` --reffa {input.fa} &> {log}
        """

