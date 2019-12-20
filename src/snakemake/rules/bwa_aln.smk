rule bwa_aln_blueprint_GRCh38:
    """
    Created:
        2018-03-16 13:02:33
    Aim:
        Align our fastq the way it is done in the Blueprint project.
    Test:
        out/bwa/aln_q-5_fa-GRCh38-Blueprint/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.sai
    Note:
        http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38
        bwa aln -q 5 grch38.fa input.fastq.gz > intermediate.sai ; bwa samse -r "read group information" grch38.fa intermediate.sai input.fastq.gz | samtools view -bS - > output.bam

    Need to process:
        inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6
        inp/fastq/run170/S001583_TH125_LC_K27ac-41049/TH125-LC-K27ac_S9
        inp/fastq/run170/S001586_TH125_SP4_K27ac-41052/TH125-SP4-K27ac_S12

        out/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.fastq
        out/gunzip/merge_illumina_fastq_sets/fastq/blueprint/8580_LC_TH101_H3K27ac/lane2_8580_ACAGTG_L002_R1.fastq
    """
    input:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["amb","ann","bwt","pac","sa"]),
        fastq="out/{filler}.fastq"
    output:
        sai="out/bwa/aln_q-{q}_{fa_genome_id}/{filler}.sai"
        #sam="out/bwa/mem_pe_{fa_genome_id}/{filler}.sam"
    params:
        prefix="out/bwa/index/{fa_genome_id}"
    wildcard_constraints:
        q="[0-9]+"
    conda:
        "../envs/bwa.yaml"
    threads:
        MAX_THREADS
    shell:
        """
        bwa\
            aln\
            -t {threads}\
            -q {wildcards.q}\
            {params.prefix}\
            {input.fastq} > {output.sai}
        """

rule bwa_samse:
    """
    Created:
        2018-03-16 13:02:33
    Aim:
        Align our fastq the way it is done in the Blueprint project.
    Test:
        out/bwa/samse_q-5_fa-GRCh38-Blueprint/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.sam
    Note:
        http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38
        bwa aln -q 5 grch38.fa input.fastq.gz > intermediate.sai ; bwa samse -r "read group information" grch38.fa intermediate.sai input.fastq.gz | samtools view -bS - > output.bam

    Need to process:
        inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6
        inp/fastq/run170/S001583_TH125_LC_K27ac-41049/TH125-LC-K27ac_S9
        inp/fastq/run170/S001586_TH125_SP4_K27ac-41052/TH125-SP4-K27ac_S12

        out/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.fastq
        out/gunzip/merge_illumina_fastq_sets/fastq/blueprint/8580_LC_TH101_H3K27ac/lane2_8580_ACAGTG_L002_R1.fastq
    """
    input:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["amb","ann","bwt","pac","sa"]),
        fastq="out/{filler}.fastq",
        sai="out/bwa/aln_q-{q}_{fa_genome_id}/{filler}.sai"
    output:
        sam="out/bwa/samse_q-{q}_{fa_genome_id}/{filler}.sam"
    params:
        prefix="out/bwa/index/{fa_genome_id}"
    conda:
        "../envs/bwa.yaml"
    wildcard_constraints:
        q="[0-9]+"
    shell:
        """
        bwa\
            samse\
            -r "@RG\\tID:IHaveNoIdea\\tSM:ofWhatToPutHere"\
            {params.prefix}\
            {input.sai}\
            {input.fastq} > {output.sam}
        """


