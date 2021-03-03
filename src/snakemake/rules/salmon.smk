rule salmon_index:
    """
    Test:
        out/salmon/index/wget/http/ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all
    """
    input:
        "out/{filler}.fa.gz"
    output:
        directory("out/salmon/index/{filler}")
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i `dirname {output}`"

rule salmon_quant_bam:
    """
    Warning: Uncomplete rule. Work todo
    Test:
        out/salmon/index/wget/http/ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all
    """
    input:
        "out/{filler}.fa.gz"
    output:
        directory("out/salmon/quant/{filler}")
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i {output}"

rule salmon_quant_fastq_pe:
    """
    Test:
        out/salmon/quant_fastq_pe_-l_ISF_--validateMappings_salmon-index-GRCh38-ensembl-r103-cdna/bedtools/bamtofastq_pe/agent/locatit_mbc_-i_-R/picard/SortSam_sortOrder-queryname/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60
    """
    input:
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz"
    output:
        directory("out/{tool}{extra}_{index_id}/{filler}")
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "salmon/quant_fastq_pe",
        index_id = "salmon-[a-zA-Z0-9-]+"
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -i {input.index} {params.extra} -1 {input.fq1} -2 {input.fq2} -o {output}"
