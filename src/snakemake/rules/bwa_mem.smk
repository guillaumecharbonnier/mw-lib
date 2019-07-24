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
        16
    shell:
        """
        bwa \
            mem \
            -t {threads} \
            {params.idxbase} \
            {input.mate1} \
            {input.mate2} > {output.sam}
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

