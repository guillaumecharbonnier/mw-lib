rule spades_se:
    """
    Created:
        2020-01-20 09:52:58
    Aim:
        Try an alternative to Edena when it fails to run because the input file is to big, e.g. T11C_all_DNA_samples without filtering.
    Doc:
        http://cab.spbu.ru/files/release3.14.0/manual.html
    Test:
        out/spades/se/ln/alias/sst/all_samples/fastq/T11C_H3K27ac/done
    """
    input:
        fq="out/{filler}.fastq.gz"
    output:
        expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["contigs.fasta", "scaffolds.fasta"]) #Many more files here
    params:
        outdir="out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="spades/se"
    conda:
        "../envs/spades.yaml"
    threads:
        16
    shell:
        """
        spades.py -t {threads} {params.extra} -s {input.fq} -o {params.outdir}
        """

rule meta_all_rnaspades_tests:
    """
    Created:
        2020-01-30 02:25:48
    """
    input:
        expand("out/rnaspades/pe_test{test}/transcripts.fasta", test=[
            "1_c_h2al2_ko_rep1", 
            "2_c_h2al2_ko",
            "3_c_h2al2",
            "4_h2al2",
            "5_h2al2/sickle/pe_-t_sanger_-q_20"])


rule rnaspades_pe_test1_c_h2al2_ko_rep1:
    """
    Created:
        2020-01-30 01:12:54
    Doc:
        http://cab.spbu.ru/files/release3.14.0/rnaspades_manual.html
    Test:
        out/rnaspades/pe_test1_c_h2al2_ko_rep1/transcripts.fasta
    """
    input:
        c_ko_r1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R1.fastq.gz",
        c_ko_r2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R2.fastq.gz"
    output:
        "out/{tool}{extra}/transcripts.fasta"
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    params:
        outdir = "out/{tool}{extra}",
        extra = params_extra
    wildcard_constraints:
        tool="rnaspades/pe_test1_c_h2al2_ko_rep1"
    conda:
        "../envs/spades.yaml"
    shell:
        """
        rnaspades.py -t {threads} {params.extra} \
            --pe1-1 {input.c_ko_r1} --pe1-2 {input.c_ko_r2} \
            -o {params.outdir}
        """


rule rnaspades_pe_test2_c_h2al2_ko:
    """
    Created:
        2020-01-30 01:12:54
    Doc:
        http://cab.spbu.ru/files/release3.14.0/rnaspades_manual.html
    Test:
        out/rnaspades/pe_test2_c_h2al2_ko/transcripts.fasta
    """
    input:
        pe1_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R1.fastq.gz",
        pe1_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R2.fastq.gz",
        pe2_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_R1.fastq.gz",
        pe2_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_R2.fastq.gz",
        pe3_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_R1.fastq.gz",
        pe3_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_R2.fastq.gz"
    output:
        "out/{tool}{extra}/transcripts.fasta"
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    params:
        outdir = "out/{tool}{extra}",
        extra = params_extra
    wildcard_constraints:
        tool="rnaspades/pe_test2_c_h2al2_ko"
    conda:
        "../envs/spades.yaml"
    shell:
        """
        rnaspades.py -t {threads} {params.extra} \
            --pe1-1 {input.pe1_1} --pe1-2 {input.pe1_2} \
            --pe2-1 {input.pe2_1} --pe2-2 {input.pe2_2} \
            --pe3-1 {input.pe3_1} --pe3-2 {input.pe3_2} \
            -o {params.outdir}
        """

rule rnaspades_pe_test3_c_h2al2:
    """
    Created:
        2020-01-30 01:12:54
    Doc:
        http://cab.spbu.ru/files/release3.14.0/rnaspades_manual.html
    Test:
        out/rnaspades/pe_test3_c_h2al2/transcripts.fasta
    """
    input:
        pe1_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R1.fastq.gz",
        pe1_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R2.fastq.gz",
        pe2_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_R1.fastq.gz",
        pe2_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_R2.fastq.gz",
        pe3_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_R1.fastq.gz",
        pe3_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_R2.fastq.gz",
        pe4_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep1_R1.fastq.gz",
        pe4_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep1_R2.fastq.gz",
        pe5_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep2_R1.fastq.gz",
        pe5_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep2_R2.fastq.gz",
        pe6_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep3_R1.fastq.gz",
        pe6_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep3_R2.fastq.gz"
    output:
        "out/{tool}{extra}/transcripts.fasta"
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    params:
        outdir = "out/{tool}{extra}",
        extra = params_extra
    wildcard_constraints:
        tool="rnaspades/pe_test3_c_h2al2"
    conda:
        "../envs/spades.yaml"
    shell:
        """
        rnaspades.py -t {threads} {params.extra} \
            --pe1-1 {input.pe1_1} --pe1-2 {input.pe1_2} \
            --pe2-1 {input.pe2_1} --pe2-2 {input.pe2_2} \
            --pe3-1 {input.pe3_1} --pe3-2 {input.pe3_2} \
            --pe4-1 {input.pe4_1} --pe4-2 {input.pe4_2} \
            --pe5-1 {input.pe5_1} --pe5-2 {input.pe5_2} \
            --pe6-1 {input.pe6_1} --pe6-2 {input.pe6_2} \
            -o {params.outdir}
        """

rule rnaspades_pe_test4_h2al2:
    """
    Created:
        2020-01-30 01:12:54
    Doc:
        http://cab.spbu.ru/files/release3.14.0/rnaspades_manual.html
    Warning:
        Does not work because RNASpades accepts a maximum of 9 pe libraries...
        Concatenation seems mandatory first.
    Test:
        out/rnaspades/pe_test4_h2al2/transcripts.fasta
    """
    input:
        pe1_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R1.fastq.gz",
        pe1_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_R2.fastq.gz",
        pe2_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_R1.fastq.gz",
        pe2_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_R2.fastq.gz",
        pe3_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_R1.fastq.gz",
        pe3_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_R2.fastq.gz",
        pe4_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep1_R1.fastq.gz",
        pe4_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep1_R2.fastq.gz",
        pe5_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep2_R1.fastq.gz",
        pe5_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep2_R2.fastq.gz",
        pe6_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep3_R1.fastq.gz",
        pe6_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep3_R2.fastq.gz",
        pe7_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep1_R1.fastq.gz",
        pe7_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep1_R2.fastq.gz",
        pe8_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep2_R1.fastq.gz",
        pe8_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep2_R2.fastq.gz",
        pe9_1 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep3_R1.fastq.gz",
        pe9_2 ="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep3_R2.fastq.gz",
        pe10_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep1_R1.fastq.gz",
        pe10_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep1_R2.fastq.gz",
        pe11_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep2_R1.fastq.gz",
        pe11_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep2_R2.fastq.gz",
        pe12_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep3_R1.fastq.gz",
        pe12_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep3_R2.fastq.gz",
        pe13_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep1_R1.fastq.gz",
        pe13_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep1_R2.fastq.gz",
        pe14_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep2_R1.fastq.gz",
        pe14_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep2_R2.fastq.gz",
        pe15_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep3_R1.fastq.gz",
        pe15_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep3_R2.fastq.gz",
        pe16_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep1_R1.fastq.gz",
        pe16_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep1_R2.fastq.gz",
        pe17_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep2_R1.fastq.gz",
        pe17_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep2_R2.fastq.gz",
        pe18_1="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep3_R1.fastq.gz",
        pe18_2="out/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep3_R2.fastq.gz",
    output:
        "out/{tool}{extra}/transcripts.fasta"
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    params:
        outdir = "out/{tool}{extra}",
        extra = params_extra
    wildcard_constraints:
        tool="rnaspades/pe_test4_h2al2"
    conda:
        "../envs/spades.yaml"
    shell:
        """
        rnaspades.py -t {threads} {params.extra} \
            --pe1-1 {input.pe1_1}   --pe1-2 {input.pe1_2} \
            --pe2-1 {input.pe2_1}   --pe2-2 {input.pe2_2} \
            --pe3-1 {input.pe3_1}   --pe3-2 {input.pe3_2} \
            --pe4-1 {input.pe4_1}   --pe4-2 {input.pe4_2} \
            --pe5-1 {input.pe5_1}   --pe5-2 {input.pe5_2} \
            --pe6-1 {input.pe6_1}   --pe6-2 {input.pe6_2} \
            --pe7-1 {input.pe7_1}   --pe7-2 {input.pe7_2} \
            --pe8-1 {input.pe8_1}   --pe8-2 {input.pe8_2} \
            --pe9-1 {input.pe9_1}   --pe9-2 {input.pe9_2} \
            --pe10-1 {input.pe10_1} --pe10-2 {input.pe10_2} \
            --pe11-1 {input.pe11_1} --pe11-2 {input.pe11_2} \
            --pe12-1 {input.pe12_1} --pe12-2 {input.pe12_2} \
            --pe13-1 {input.pe13_1} --pe13-2 {input.pe13_2} \
            --pe14-1 {input.pe14_1} --pe14-2 {input.pe14_2} \
            --pe15-1 {input.pe15_1} --pe15-2 {input.pe15_2} \
            --pe16-1 {input.pe16_1} --pe16-2 {input.pe16_2} \
            --pe17-1 {input.pe17_1} --pe17-2 {input.pe17_2} \
            --pe18-1 {input.pe18_1} --pe18-2 {input.pe18_2} \
            -o {params.outdir}
        """


rule rnaspades_pe_test5_h2al2_filler:
    """
    Created:
        2020-01-30 01:12:54
    Aim:
        Testing RNASpades with our without read prefiltering
    Doc:
        http://cab.spbu.ru/files/release3.14.0/rnaspades_manual.html
    Test:
        out/rnaspades/pe_test5_h2al2/sickle/pe_-t_sanger_-q_20/transcripts.fasta
    """
    input:
        pe1_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_1.fastq.gz",
        pe1_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep1_2.fastq.gz",
        pe2_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_1.fastq.gz",
        pe2_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep2_2.fastq.gz",
        pe3_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_1.fastq.gz",
        pe3_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-KO-Rep3_2.fastq.gz",
        pe4_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep1_1.fastq.gz",
        pe4_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep1_2.fastq.gz",
        pe5_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep2_1.fastq.gz",
        pe5_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep2_2.fastq.gz",
        pe6_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep3_1.fastq.gz",
        pe6_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-C-H2AL2-WT-Rep3_2.fastq.gz",
        pe7_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep1_1.fastq.gz",
        pe7_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep1_2.fastq.gz",
        pe8_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep2_1.fastq.gz",
        pe8_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep2_2.fastq.gz",
        pe9_1 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep3_1.fastq.gz",
        pe9_2 ="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-KO-Rep3_2.fastq.gz",
        pe10_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep1_1.fastq.gz",
        pe10_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep1_2.fastq.gz",
        pe11_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep2_1.fastq.gz",
        pe11_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep2_2.fastq.gz",
        pe12_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep3_1.fastq.gz",
        pe12_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-R-H2AL2-WT-Rep3_2.fastq.gz",
        pe13_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep1_1.fastq.gz",
        pe13_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep1_2.fastq.gz",
        pe14_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep2_1.fastq.gz",
        pe14_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep2_2.fastq.gz",
        pe15_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep3_1.fastq.gz",
        pe15_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-KO-Rep3_2.fastq.gz",
        pe16_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep1_1.fastq.gz",
        pe16_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep1_2.fastq.gz",
        pe17_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep2_1.fastq.gz",
        pe17_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep2_2.fastq.gz",
        pe18_1="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep3_1.fastq.gz",
        pe18_2="out/{filler}/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-P-H2AL2-WT-Rep3_2.fastq.gz"
    output:
        "out/{tool}{extra}/{filler}/transcripts.fasta"
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="rnaspades/pe_test5_h2al2"
    conda:
        "../envs/spades.yaml"
    shell:
        """
        rnaspades.py -t {threads} {params.extra} \
            --pe1-1  {input.pe1_1}  --pe1-2  {input.pe1_2} \
            --pe2-1  {input.pe2_1}  --pe2-2  {input.pe2_2} \
            --pe3-1  {input.pe3_1}  --pe3-2  {input.pe3_2} \
            --pe4-1  {input.pe4_1}  --pe4-2  {input.pe4_2} \
            --pe5-1  {input.pe5_1}  --pe5-2  {input.pe5_2} \
            --pe6-1  {input.pe6_1}  --pe6-2  {input.pe6_2} \
            --pe7-1  {input.pe7_1}  --pe7-2  {input.pe7_2} \
            --pe8-1  {input.pe8_1}  --pe8-2  {input.pe8_2} \
            --pe9-1  {input.pe9_1}  --pe9-2  {input.pe9_2} \
            --pe10-1 {input.pe10_1} --pe10-2 {input.pe10_2} \
            --pe11-1 {input.pe11_1} --pe11-2 {input.pe11_2} \
            --pe12-1 {input.pe12_1} --pe12-2 {input.pe12_2} \
            --pe13-1 {input.pe13_1} --pe13-2 {input.pe13_2} \
            --pe14-1 {input.pe14_1} --pe14-2 {input.pe14_2} \
            --pe15-1 {input.pe15_1} --pe15-2 {input.pe15_2} \
            --pe16-1 {input.pe16_1} --pe16-2 {input.pe16_2} \
            --pe17-1 {input.pe17_1} --pe17-2 {input.pe17_2} \
            --pe18-1 {input.pe18_1} --pe18-2 {input.pe18_2} \
            -o {params.outdir}
        """


rule rnaspades_1_pe_lib:
    """
    Created:
        2020-01-30 01:12:54
    Aim:
        Testing RNASpades with our without read prefiltering
    Doc:
        http://cab.spbu.ru/files/release3.14.0/rnaspades_manual.html
    Test:
        cat-h2al2-rna_1.fastq.gz    expand("sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}_1.fastq.gz", stage=["P","R","C"], condition=["WT","KO"], replicate=["1","2","3"])
        cat-h2al2-rna_2.fastq.gz    expand("sickle/pe_-t_sanger_-q_20/ln/pe_remove_mate_prefix/cat/merge_lanes_nextseq500_pe_LN_RM/ln/updir/mw-sk/inp/fastq/tgml/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}_2.fastq.gz", stage=["P","R","C"], condition=["WT","KO"], replicate=["1","2","3"])

        out/rnaspades/1_pe_lib/cat/cat-h2al2-rna/transcripts.fasta
    """
    input:
        pe1_1 ="out/{filler}_1.fastq.gz",
        pe1_2 ="out/{filler}_2.fastq.gz"
    output:
        "out/{tool}{extra}/{filler}/transcripts.fasta"
        #expand("out/{{tool}}{{extra}}/{{filler}}/{filename}",filename=["transcripts.fasta"]) # Many more files here
    benchmark:
        "out/{tool}{extra}/{filler}/transcripts.benchmark.tsv"
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="rnaspades/1_pe_lib"
    conda:
        "../envs/spades.yaml"
    threads:
        16
    shell:
        """
        rnaspades.py -t {threads} -m 100 {params.extra} \
            --pe1-1  {input.pe1_1}  --pe1-2  {input.pe1_2} \
            -o {params.outdir}
        """


