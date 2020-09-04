rule MindTheGap_find_se:
    """
    Created:
        2020-08-18 20:54:37
    Test:
        out/MindTheGap/find_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.done
    """
    input:
        fastq="out/{filler}.fastq",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        breakpoints="out/{tool}_{fa_genome_id}/{filler}.breakpoints",
        h5="out/{tool}_{fa_genome_id}/{filler}.h5",
        vcf="out/{tool}_{fa_genome_id}/{filler}.othervariants.vcf"
    params:
        outdir= "out/{tool}_{fa_genome_id}/{filler}"
    wildcard_constraints:
        tool = "MindTheGap/find_se"
    conda:
        "../envs/mindthegap.yaml"
    shell:
        """
        MindTheGap find -in {input.fastq} -ref {input.fasta} -out {params.outdir}
        """

rule MindTheGap_find_pe:
    """
    Created:
        2020-08-18 20:54:37
    Test:
        out/MindTheGap/fill/MindTheGap/find_pe_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/856_H3K27ac.insertions.vcf
    """
    input:
        fastq_1="out/{filler}_1.fastq",
        fastq_2="out/{filler}_2.fastq",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        breakpoints="out/{tool}_{fa_genome_id}/{filler}.breakpoints",
        h5="out/{tool}_{fa_genome_id}/{filler}.h5",
        vcf="out/{tool}_{fa_genome_id}/{filler}.othervariants.vcf"
    params:
        outdir= "out/{tool}_{fa_genome_id}/{filler}"
    wildcard_constraints:
        tool = "MindTheGap/find_pe"
    conda:
        "../envs/mindthegap.yaml"
    shell:
        """
        MindTheGap find -in {input.fastq_1},{input.fastq_2} -ref {input.fasta} -out {params.outdir}
        """

#fill build/bin/MindTheGap fill -graph example.h5 -bkpt example.breakpoints -out example
# 3 files are generated:
#   example.insertions.fasta (insertion sequences)
#   example.insertions.vcf (insertion variants)
#   example.info.txt (log file)

rule MindTheGap_fill:
    """
    
    Test:
        out/MindTheGap/fill/MindTheGap/find_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.insertions.vcf
    """
    input:
        h5   = "out/{filler}.h5",
        bkpt = "out/{filler}.breakpoints"
    output:
        insertions_fasta = "out/{tool}/{filler}.insertions.fasta",
        insertions_vcf   = "out/{tool}/{filler}.insertions.vcf"
    params:
        outdir = "out/{tool}/{filler}"
    wildcard_constraints:
        tool = "MindTheGap/fill"
    conda:
        "../envs/mindthegap.yaml"
    shell:
        """
        MindTheGap fill -graph {input.h5} -bkpt {input.bkpt} -out {params.outdir}
        """


