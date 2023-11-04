rule salmon_quant_bam:
    """
    Warning: Uncomplete rule. Work todo
    Test:
        out/salmon/index/wget/http/ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all
    """
    input:
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        bam = "out/{filler}.bam"
    output:
        directory("out/{tool}{extra}_{index_id}/{filler}")
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "salmon/quant_bam",
        index_id = "salmon-[a-zA-Z0-9-]+"
    threads:
        MAX_THREADS
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -p {threads} -i {input.index} {params.extra} -a {input.bam} -o {output}"

rule salmon_quant_fastq_pe:
    """
    Test:
        out/salmon/quant_fastq_pe_-l_A_--validateMappings_salmon-index-GRCh38-ensembl-r102/sickle/pe_-t_sanger_-q_20/sra-tools/fastq-dump_pe/SRR4123954/quant.sf
    """
    input:
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz"
    output:
        "out/{tool}{extra}_{index_id}/{filler}/quant.sf"
        #directory("out/{tool}{extra}_{index_id}/{filler}")
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "salmon/quant_fastq_pe",
        index_id = "salmon-[a-zA-Z0-9-]+"
    threads:
        MAX_THREADS
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -p {threads} -i {input.index} {params.extra} -1 {input.fq1} -2 {input.fq2} -o `dirname {output}`"

rule salmon_quant_fastq_se:
    """
    Test:
        out/salmon/quant_fastq_se_-l_A_--validateMappings_--numGibbsSamples_20_--gcBias_salmon-index-GRCh38-ensembl-r102/sickle/se_-t_sanger_-q_20/sra-tools/fastq-dump_se/SRR1927116/quant.sf
    """
    input:
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        fq = "out/{filler}.fastq.gz",
    output:
        "out/{tool}{extra}_{index_id}/{filler}/quant.sf"
        #directory("out/{tool}{extra}_{index_id}/{filler}")
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "salmon/quant_fastq_se",
        index_id = "salmon-[a-zA-Z0-9-]+"
    threads:
        MAX_THREADS
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -p {threads} -i {input.index} {params.extra} -r {input.fq} -o `dirname {output}`"

