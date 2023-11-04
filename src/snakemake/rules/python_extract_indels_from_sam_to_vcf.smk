rule python_extract_indels_from_sam_to_vcf:
    """
    Created:
        2020-08-22 15:01:15
    Test:
        out/python/extract_indels_from_sam_to_vcf_fa-genome-GRCh38-ensembl-r100/bowtie2/se_fa_Abraham2017_GRCh38-ensembl-r100/megahit/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.contigs.vcf
    """
    input:
        sam = "out/{filler}.sam",
        fa = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id]),
        py = "../mw-lib/src/python/extract_indels_from_sam_to_vcf.py"
    output:
        vcf = "out/python/extract_indels_from_sam_to_vcf_{fa_genome_id}/{filler}.vcf"
    log:
        "out/python/extract_indels_from_sam_to_vcf_{fa_genome_id}/{filler}.log"
    conda:
        "../envs/pysam.yaml"
    shell:
        """
        {input.py} -s {input.sam} -g {input.fa} -v {output.vcf} 2> {log}
        """
