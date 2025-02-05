rule phantompeakqualtools_bam_noctrl:
    """
    Created:
        2019-02-12 10:49:13
    Doc:
        https://code.google.com/archive/p/phantompeakqualtools/
    Warning:
        It is EXTREMELY important to filter out multi-mapping reads from the BAM/tagAlign files. Large number of multimapping reads can severly affect the phantom peak coefficient and peak calling results.
        If a dataset seems to have high PCR bottlenecking, then you might want to actually clamp the number of unique mappping reads per position to 1 or upto 5. If not the phantom peak coefficient can be artificially good.
    Note:
        multiqc expect '*spp.out'
    Test:
        out/phantompeakqualtools/bam_noctrl_-savp/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243.tsv
    """
    input:
        bam = "out/{filler}.bam"
    output:
        tsv     = "out/{tool}{extra}/{filler}.tsv",
        multiqc = "out/{tool}{extra}/{filler}.spp.out",
        pdf     = "out/{tool}{extra}/{filler}.pdf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="phantompeakqualtools/bam_noctrl"
    conda:
        "../envs/phantompeakqualtools.yaml"
    shell:
        "run_spp.R -c={input.bam} {params} -odir=`dirname {output.tsv}` -out={output.tsv}; "
        "ln -srf {output.tsv} {output.multiqc}"

rule phantompeakqualtools_bam_ctrl:
    """
    Created:
        2019-02-12 10:49:13
    Doc:
        https://code.google.com/archive/p/phantompeakqualtools/
    Warning:
        It is EXTREMELY important to filter out multi-mapping reads from the BAM/tagAlign files. Large number of multimapping reads can severly affect the phantom peak coefficient and peak calling results.
        If a dataset seems to have high PCR bottlenecking, then you might want to actually clamp the number of unique mappping reads per position to 1 or upto 5. If not the phantom peak coefficient can be artificially good.
    Test:
        out/phantompeakqualtools/bam_ctrl_idr/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_VS_SRR3126242.txt
    """
    input:
        bam = "out/{filler}/{bam}.bam",
        ctrl= "out/{filler}/{ctrl}.bam"
    output:
        #"out/{tool}{extra}/{filler}.done"
        "out/{tool}{extra}/{filler}/{bam}_VS_{ctrl}.txt"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="phantompeakqualtools/bam_ctrl"
    conda:
        "../envs/phantompeakqualtools.yaml"
    shell:
        "run_spp.R -c={input.bam} -i={input.ctrl} {params} -odir=`dirname {output}` -out={output}"

