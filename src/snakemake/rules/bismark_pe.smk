rule bismark_tests:
    """
    Note about the target samples:
        - 728 because it is the smallest file
        - 525 because it is the largest file that may see reads for spike-in
    """
    input:
        # Note the -q_15 added to trim_galore in an attempt to keep R2 intact for downstream and maybe improve the mapping efficiency.
        expand(
            "out/bismark/pe_--pbat_{genome}/trim_galore/pe_-q_15_--non_directional_--rrbs_--length_15_--fastqc/ln/pe_remove_mate_prefix/ln/updir/mw-tall-rrbs/inp/bcl/230725_NB501645_0730_AHMTMKBGXT/Data/Intensities/BaseCalls/fumi_tools_demultiplex/{sample}.bam",
            genome = ["fa-genome-GRCh38-main-chr", "fa-genome-RRBS-v2-methylated-control", "fa-genome-RRBS-v2-unmethylated-control"],
            sample = ["Undetermined_S0_L001", "728-BARLYD_S728-BARLYD_L001", "525-HALJEA_S525-HALJEA_L001"]
        ),
        expand(
            "out/bismark/pe_--pbat_{genome}/trim_galore/pe_--non_directional_--rrbs_--length_15_--fastqc/ln/updir/mw-tall-rrbs/inp/diagenode_test_fastq/RRBS_V2_sample_test.bam",
            genome = ["fa-genome-GRCh38-main-chr", "fa-genome-RRBS-v2-methylated-control", "fa-genome-RRBS-v2-unmethylated-control"]
        )

rule bismark_pe:
    """
    Run bismark on paired-end reads.

    Example for Diagenode RRBS V2:
        # bismark -q --pbat --prefix Meth_ctrl ./genomes/RRBS_methylated_control -1 MySample_R1_val_1.fastq -2 MySample_R2_val_2.fastq 
    """
    input:
        fwd = "out/{filler}_1.fastq.gz",
        rev = "out/{filler}_2.fastq.gz",
        # ref = "out/bismark/genome_preparation/{fa_genome_id}/Bisulfite_Genome"
        ref = directory("out/bismark/genome_preparation/{fa_genome_id}/Bisulfite_Genome")
    output:
        bam = "out/{tool}{extra}_{fa_genome_id}/{filler}.bam"
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}.log"
    params:
        extra = params_extra,
        tmp_dir = "out/{tool}{extra}_{fa_genome_id}/{filler}_tmp"
    conda:
        "../envs/bismark.yaml"
    wildcard_constraints:
        tool = "bismark/pe"
    resources:
        high_io = 1,
        ram = 40 # This is without using the --multicore option for human genome.
    threads: 
        2 # Bismark will run multiple tasks in parallel even without settings --multicore. See https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
    shell:
        """
        bismark {params.extra} \
            -o $(dirname {output.bam}) $(dirname {input.ref}) \
            --temp_dir {params.tmp_dir} \
            -1 {input.fwd} -2 {input.rev} &> {log}
        mv out/{wildcards.tool}{wildcards.extra}_{wildcards.fa_genome_id}/{wildcards.filler}_1_bismark_bt2_pe.bam {output.bam}
        """
