rule bedtools_genomecov_bam_to_d:
    """
    Created:
        2019-09-05 16:08:03
    Doc:
        https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html#d-reporting-per-base-genome-coverage
    
    """
    input:
        bam="out/{filler}.bam"
    output:
        d="out/{tool}{extra}/{filler}.d"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/genomecov_bam_to_d"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools genomecov {params} -d -ibam {input.bam} > {output.d}"


rule bedtools_genomecov_bam_to_bga:
    """
    Created:
        2018-03-06 18:03:38

    Note:
        2. If the input is in BAM (-ibam) format, the BAM file must be sorted by position. Using samtools sort aln.bam aln.sorted will suffice.
        *****WARNING: Genome (-g) files are ignored when BAM input is provided. 
    Test:
        out/bedtools/genomecov_bam_to_bga/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bg
    """
    input:
        bam="out/{filler}.bam"
    output:
        bedgraph="out/{tool}{extra}/{filler}.bg"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/genomecov_bam_to_bga"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools genomecov {params} -bga -ibam {input.bam} > {output.bedgraph}"

rule bedtools_genomecov_bed_to_bga:
    """
    Created:
        2018-03-06 18:03:38
    Note:
         If using BED/GFF/VCF, the input (-i) file must be grouped by chromosome. A simple sort -k 1,1 in.bed > in.sorted.bed will suffice. Also, if using BED/GFF/VCF, one must provide a genome file via the -g argument.
    Note:
        Using "bedgraph" as extension because it is the one recognized by IGV:
        https://software.broadinstitute.org/software/igv/bedgraph
    Test:
        out/bedtools/genomecov_bed_to_bga_chrominfo-mm10/bedtools/bamtobed/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bg
    """
    input:
        bam = "out/{filler}.bed",
        chromInfo = lambda wildcards: eval(config['ids'][wildcards.chrominfo_id])
    output:
        bedgraph="out/{tool}{extra}_{chrominfo_id}/{filler}.bg"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bedtools/genomecov_bed_to_bga"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools genomecov {params} -bga -i {input.bam} -g {input.chromInfo} > {output.bedgraph}"


