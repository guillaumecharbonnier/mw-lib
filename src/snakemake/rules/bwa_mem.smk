rule bwa_mem_pe:
    """
    Created:
        2017-05-14 13:14:45
    Test:
        out/bwa/mem_pe_GRCm38/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.sam
    """
    input:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["amb","ann","bwt","pac","sa"]),
        mate1="out/{filler}_1.fastq",
        mate2="out/{filler}_2.fastq"
    output:
        sam="out/bwa/mem_pe_{fa_genome_id}/{filler}.sam"
    params:
        idxbase="out/bwa/index/{fa_genome_id}"
    conda:
        "../envs/bwa.yaml"
    threads:
        MAX_THREADS
    shell:
        """
        bwa \
            mem \
            -t {threads} \
            {params.idxbase} \
            {input.mate1} \
            {input.mate2} > {output.sam}
        """

rule bwa_mem_se:
    """
    Aim:
        Currently rewriting this rule in order to be able to ask for bowtie index with alt haplotype on demand.
    Created:
        2019-12-01 23:22:39
    Test:
        out/bwa/mem2_se_bwa-index-hg19-main-chr-and-contigs-with-inserts-from-T11C-H3K27ac/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac.sam
    """
    input:
        fastq="out/{filler}.fastq",
        index_files=lambda wildcards: eval(config['ids'][wildcards.bwa_index_id])
    output:
        sam="out/bwa/mem2_se_{bwa_index_id}/{filler}.sam"
    #params:
        #extra=params_extra,
        #idxbase=params_bwa_index_base_path
    conda:
        "../envs/bwa.yaml"
    threads:
        MAX_THREADS
    shell:
        """
        the_first_index_file=`echo {input.index_files} | cut -f1 -d' ' `

        echo $the_first_index_file

        the_first_index_file_without_ext=${{the_first_index_file%.*}}

        echo $the_first_index_file_without_ext

        bwa mem -t {threads} $the_first_index_file_without_ext {input.fastq} > {output.sam}
        """

rule bwa_mem_se_legacy:
    """
    Created:
        2019-12-01 23:22:39
    Test:
        out/bwa/mem_se_fa-genome-GRCh38/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/Jurkat_SRR1057274_H3K27ac.sam
    """
    input:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["amb","ann","bwt","pac","sa"]),
        fastq="out/{filler}.fastq",
    output:
        sam="out/bwa/mem_se_{fa_genome_id}/{filler}.sam"
    params:
        idxbase="out/bwa/index/{fa_genome_id}"
    conda:
        "../envs/bwa.yaml"
    threads:
        MAX_THREADS
    shell:
        """
        bwa \
            mem \
            -t {threads} \
            {params.idxbase} \
            {input.fastq} > {output.sam}
        """


"""
        sam="out/bowtie2/pe_{index}/{filler}.sam",
        unmapped_single="out/bowtie2/pe_{index}/{filler}/unmapped/single.fastq.gz",
        unmapped_pair1="out/bowtie2/pe_{index}/{filler}/unmapped/pair.fastq.1.gz",
        unmapped_pair2="out/bowtie2/pe_{index}/{filler}/unmapped/pair.fastq.2.gz",
        # alias because the default extension for unmapped is not convenient for 'fastq' suffix handling.
        ln_unmapped_single="out/bowtie2/pe_{index}/{filler}/unmapped_single.fastq.gz",
        ln_unmapped_pair1="out/bowtie2/pe_{index}/{filler}/unmapped_pair1.fastq.gz",
        ln_unmapped_pair2="out/bowtie2/pe_{index}/{filler}/unmapped_pair2.fastq.gz",
        bowtie2log="out/bowtie2/pe_{index}/{filler}.log"
    params:
        unmapped_pair_base="out/bowtie2/pe_{index}/{filler}/unmapped/pair.fastq.gz",
        index=bowtie2_index_base_path
"""

