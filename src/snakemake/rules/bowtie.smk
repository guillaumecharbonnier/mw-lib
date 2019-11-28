rule bowtie_single_end_extra:
    """
    Modified:
        2017-05-06 14:53:46 - Tool is now installed with conda.
    Doc:
        http://bowtie-bio.sourceforge.net/manual.shtml
    Test:
        out/bowtie/se_--chunkmbs_256_--best_--strata_-m_1_-n_2_ebwt-hg19/sra-tools/fastq-dump_se/SRR1202037.sam
    TODO: CHECK GZ AS INPUT FOR THIS RULE?
    """
    input:
        fastq="out/{filler}.fastq.gz",
        index_files= lambda wildcards: eval(config['ids'][wildcards.ebwt_id])
    output:
        sam      = "out/{tool}{extra}_{ebwt_id}/{filler}.sam",
        unmapped = "out/{tool}{extra}_{ebwt_id}/{filler}_unmapped.fastq.gz",
    log:
                   "out/{tool}{extra}_{ebwt_id}/{filler}.log"
    benchmark:
                   "out/{tool}{extra}_{ebwt_id}/{filler}.benchmark.tsv"
    params:
        extra=params_extra,
        #index_basepath=params_bowtie2_index_base_path
        #index_basepath="hg19" 
        index_basepath= lambda wildcards: os.path.commonprefix(eval(config['ids'][wildcards.ebwt_id]))[:-1]
    wildcard_constraints:
        tool="bowtie/se"
    conda:
        "../envs/bowtie.yaml"
    threads:
        16
    shell:
        """
        ( UN=`echo {output.unmapped} | sed 's/.gz$//'`
        gzip -dc {input.fastq} | bowtie {params.extra} -p {threads} -S  {params.index_basepath} --un $UN --max /dev/null - > {output.sam}
        gzip $UN ) &> {log}
        """
        #bowtie {params.extra} -p {threads} -S  {params.index_basepath} --un  {output.unmapped_single} --max /dev/null {input.fastq} > {output.sam}"
        #"bowtie --chunkmbs 256 --best --strata -m 1 -n 2 -p {threads} -S  {params.index_basepath} --un  {output.unmapped_single} --max /dev/null {input.fastq} > {output.sam}"


#rule bowtie2_paired_end_extra:
#    """
#    Created:
#        2018-08-23 14:40:53
#    Aim:
#        Produce aligned read in bam format from fastq files.
#        Preset added to test for "--very-sensitive-local" in Fasktd1 project.
#    Test:
#        out/bowtie2/pe_hg19/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR3938832.sam
#    """
#    input:
#        mate1="out/{filler}_1.fastq.gz",
#        mate2="out/{filler}_2.fastq.gz",
#        index_part=input_bowtie2_index_parts
#    output:
#        sam                = temp("out/{tool}{extra}_{index}/{filler}.sam"),
#        unmapped_single    = "out/{tool}{extra}_{index}/{filler}/unmapped/single.fastq.gz",
#        unmapped_pair1     = "out/{tool}{extra}_{index}/{filler}/unmapped/pair.fastq.1.gz",
#        unmapped_pair2     = "out/{tool}{extra}_{index}/{filler}/unmapped/pair.fastq.2.gz",
#        # alias because the defaul{tool}{extra}_{index}t extension for un_{extra}mapped is not convenient for 'fastq' suffix handling.
#        ln_unmapped_single = "out/{tool}{extra}_{index}/{filler}/unmapped_single.fastq.gz",
#        ln_unmapped_pair1  = "out/{tool}{extra}_{index}/{filler}/unmapped_pair1.fastq.gz",
#        ln_unmapped_pair2  = "out/{tool}{extra}_{index}/{filler}/unmapped_pair2.fastq.gz"
#    log:
#        "out/{tool}{extra}_{index}/{filler}.log"
#    params:
#        extra = params_extra,
#        unmapped_pair_base = "out/{tool}{extra}_{index}/{filler}/unmapped/pair.fastq.gz",
#        index=params_bowtie2_index_base_path
#    conda:
#        "../envs/bowtie2.yaml"
#    wildcard_constraints:
#        tool="bowtie2/pe"
#    threads:
#        16
#    shell:
#        """
#        bowtie2 -p {threads} -x {params.index}\
#            -1 {input.mate1} -2 {input.mate2}\
#            --un-gz {output.unmapped_single}\
#            --un-conc-gz {params.unmapped_pair_base}\
#            -S {output.sam}\
#            {params.extra}\
#            2> {log}
#
#        ln -srf {output.unmapped_single} {output.ln_unmapped_single}
#        ln -srf {output.unmapped_pair1} {output.ln_unmapped_pair1}
#        ln -srf {output.unmapped_pair2} {output.ln_unmapped_pair2}
#        """
#
