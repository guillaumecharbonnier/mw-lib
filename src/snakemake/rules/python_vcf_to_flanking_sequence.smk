rule python_vcf_to_flanking_sequence:
    """
    Created:
        2020-08-22 15:01:15
    Test:
        out/python/vcf_to_flanking_sequence_fa-genome-GRCh38-ensembl-r100/grep/extract-indel-from-vcf/bcftools/mpileup_fa-genome-GRCh38-ensembl-r100/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_GRCh38-ensembl-r100/spades/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/contigs.fasta
    """
    input:
        vcf = "out/{filler}.vcf",
        fa = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id]),
        py = "../mw-lib/src/python/vcf_to_flanking_sequence.py"
    output:
        fa = "out/python/vcf_to_flanking_sequence_{fa_genome_id}/{filler}.fasta"
    shell:
        """
        {input.py} -v {input.vcf} -g {input.fa} > {output.fa}
        """
