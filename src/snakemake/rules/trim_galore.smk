rule trim_galore_pe:
    """
    Test:
        out/trim_galore/pe/ln/pe_remove_mate_prefix/ln/updir/mw-tall-rrbs/inp/bcl/230725_NB501645_0730_AHMTMKBGXT/Data/Intensities/BaseCalls/fumi_tools_demultiplex/Undetermined_S0_L001_1.fastq.gz
    """
    input:
        ["out/{filler}_1.fastq.gz", "out/{filler}_2.fastq.gz"],
    output:
        fasta_fwd  = "out/{tool}{extra}/{filler}_1.fastq.gz",
        report_fwd = "out/{tool}{extra}/{filler}_1.fastq.gz_trimming_report.txt",
        fasta_rev  = "out/{tool}{extra}/{filler}_2.fastq.gz",
        report_rev = "out/{tool}{extra}/{filler}_2.fastq.gz_trimming_report.txt",
    threads: 1
    params:
        extra = params_extra,
        # extra = "--illumina -q 20",
    log:
        "out/{tool}{extra}/{filler}.log",
    wildcard_constraints:
        tool = "trim_galore/pe"
    # Using the wrapper fails
    # wrapper:
    #     "v2.3.2/bio/trim_galore/pe"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        """
        trim_galore {params.extra} --paired --output_dir $(dirname {output.fasta_fwd}) {input[0]} {input[1]} > {log}
        # Renaming trim_galore output fastq files to match the output of the other tools
        mv out/{wildcards.tool}{wildcards.extra}/{wildcards.filler}_1_val_1.fq.gz {output.fasta_fwd}
        mv out/{wildcards.tool}{wildcards.extra}/{wildcards.filler}_2_val_2.fq.gz {output.fasta_rev}
        """