rule bowtie2_single_end_extra_misc_inps:
    """
    Modified:
        2017-05-06 14:53:46 - Tool is now installed with conda.
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
        out/bowtie2/se_mm10/sra-tools/fastq-dump/SRR1202037.sam
        out/bowtie2/se_vsl_mm10/sra-tools/fastq-dump/SRR1202037.sam
        out/bowtie2/se_vsl_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump/SRR1202037.sam
    """
    input:
        fastq="out/{filler}.{ext}",
        index_dir=input_bowtie2_index_parts
    output:
        sam             = "out/{tool}_{ext}{extra}_{index}/{filler}.sam",
        unmapped_single = "out/{tool}_{ext}{extra}_{index}/{filler}_single.fastq.gz",
    log:
                          "out/{tool}_{ext}{extra}_{index}/{filler}.log"
    params:
        extra=params_extra,
        index=params_bowtie2_index_base_path
    wildcard_constraints:
        tool="bowtie2/se",
        ext="fastq|fasta|fastq.gz|fasta.gz"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        MAX_THREADS
    shell: 
        """
        bowtie2 -p {threads} -x {params.index} \
            -U {input.fastq} \
            --un-gz {output.unmapped_single} \
            -S {output.sam} \
            {params.extra}\
            2> {log}
        """

rule bowtie2_single_end_extra:
    """
    Modified:
        2017-05-06 14:53:46 - Tool is now installed with conda.
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
        out/bowtie2/se_mm10/sra-tools/fastq-dump/SRR1202037.sam
        out/bowtie2/se_vsl_mm10/sra-tools/fastq-dump/SRR1202037.sam
        out/bowtie2/se_vsl_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump/SRR1202037.sam
    """
    input:
        fastq="out/{filler}.fastq.gz",
        index_dir=input_bowtie2_index_parts
    output:
        sam             = "out/{tool}{extra}_{index}/{filler}.sam",
        unmapped_single = "out/{tool}{extra}_{index}/{filler}_single.fastq.gz",
    log:
                          "out/{tool}{extra}_{index}/{filler}.log"
    params:
        extra=params_extra,
        index=params_bowtie2_index_base_path
    wildcard_constraints:
        tool="bowtie2/se"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        MAX_THREADS
    shell: 
        """
        bowtie2 -p {threads} -x {params.index} \
            -U {input.fastq} \
            --un-gz {output.unmapped_single} \
            -S {output.sam} \
            {params.extra}\
            2> {log}
        """

rule bowtie2_paired_end_extra:
    """
    Created:
        2018-08-23 14:40:53
    Aim:
        Produce aligned read in bam format from fastq files.
        Preset added to test for "--very-sensitive-local" in Fasktd1 project.
    Test:
        # SRR3938832 is RNA-Seq while SRR1509753 is ChIP-Seq
        #out/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR3938832.sam
        out/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR1509753.sam
    """
    input:
        mate1="out/{filler}_1.fastq.gz",
        mate2="out/{filler}_2.fastq.gz",
        index_part=input_bowtie2_index_parts
    output:
        sam                = temp("out/{tool}{extra}_{index}/{filler}.sam"),
        unmapped_single    = "out/{tool}{extra}_{index}/{filler}/unmapped/single.fastq.gz",
        unmapped_pair1     = "out/{tool}{extra}_{index}/{filler}/unmapped/pair.fastq.1.gz",
        unmapped_pair2     = "out/{tool}{extra}_{index}/{filler}/unmapped/pair.fastq.2.gz",
        # alias because the defaul{tool}{extra}_{index}t extension for un_{extra}mapped is not convenient for 'fastq' suffix handling.
        ln_unmapped_single = "out/{tool}{extra}_{index}/{filler}/unmapped_single.fastq.gz",
        ln_unmapped_pair1  = "out/{tool}{extra}_{index}/{filler}/unmapped_pair1.fastq.gz",
        ln_unmapped_pair2  = "out/{tool}{extra}_{index}/{filler}/unmapped_pair2.fastq.gz"
    log:
        "out/{tool}{extra}_{index}/{filler}.log"
    params:
        extra = params_extra,
        unmapped_pair_base = "out/{tool}{extra}_{index}/{filler}/unmapped/pair.fastq.gz",
        index=params_bowtie2_index_base_path
    conda:
        "../envs/bowtie2.yaml"
    wildcard_constraints:
        tool="bowtie2/pe"
    threads:
        MAX_THREADS
    shell:
        """
        bowtie2 -p {threads} -x {params.index}\
            -1 {input.mate1} -2 {input.mate2}\
            --un-gz {output.unmapped_single}\
            --un-conc-gz {params.unmapped_pair_base}\
            -S {output.sam}\
            {params.extra}\
            2> {log}

        ln -srf {output.unmapped_single} {output.ln_unmapped_single}
        ln -srf {output.unmapped_pair1} {output.ln_unmapped_pair1}
        ln -srf {output.unmapped_pair2} {output.ln_unmapped_pair2}
        """

