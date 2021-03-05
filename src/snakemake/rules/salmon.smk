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
        "salmon index -t {input} -i {output}"


rule salmon_index_ensembl_with_decoys:
    """
    Test:
        out/salmon/index_ensembl_with_decoys/release-103/fasta/homo_sapiens/Homo_sapiens.GRCh38
    """
    input:
        cdna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/cdna/{specie_assembly}.cdna.all.fa.gz",
        ncrna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/ncrna/{specie_assembly}.ncrna.fa.gz",
        dna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/dna/{specie_assembly}.dna.primary_assembly.fa.gz"
    output:
        directory("out/salmon/index_ensembl_with_decoys/{release_fasta_specie}/{specie_assembly}")
    wildcard_constraints:
        release_fasta_specie="release-[0-9]+/fasta/[a-z_]+",
        specie_assembly="[A-Za-z0-9_.]+"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        mkdir -p {output}
        zcat {input} > {output}/transcriptome_genome.fa

        zcat {input.dna} | grep "^>" | cut -d " " -f 1 | sed "s/>//g" > {output}/decoys.txt

        salmon index -t {output}/transcriptome_genome.fa -d {output}/decoys.txt -i {output}

        rm -f {output}/transcriptome_genome.fa.gz
        """

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
        4
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -p {threads} -i {input.index} {params.extra} -a {input.bam} -o {output}"

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
    threads:
        4
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon quant -p {threads} -i {input.index} {params.extra} -1 {input.fq1} -2 {input.fq2} -o {output}"
