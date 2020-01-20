rule bwa_index:
    """
    Created:
        2017-05-14 13:14:45
    Doc:
        http://bio-bwa.sourceforge.net/bwa.shtml
    Test:
        out/bwa/index/fa-genome-GRCh38-Blueprint.amb
    """
    input:
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
        #fai   = lambda wildcards: [path + '.fai' for path in eval(config['ids'][wildcards.bam_list_id])],
        #genome=input_fa_genome
    output:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["amb","ann","bwt","pac","sa"])
    params:
        prefix="out/bwa/index/{fa_genome_id}"
    conda:
        "../envs/bwa.yaml"
    threads:
        1
    shell:
        "bwa index -p {params.prefix} {input.fasta}"

rule bwa_index_alt:
    """
    Created:
        2020-01-17 19:10:25
    Aim:
        Create a bwa index with alt contigs

    TODO: Wait for test below to complete. If alt is not integrated, try to just add it in output folder with 'ln' before running bwa index
    Doc:
        http://bio-bwa.sourceforge.net/bwa.shtml
    Test:
        out/bwa/index_alt/awk/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_fa-genome-hg19-main-chr/bowtie2/se_fasta_Abraham2017_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_all_DNA_samples_contigs.alt
    """
    input:
        fasta = "out/{filler}.fa",
        alt   = "out/{filler}.alt"
    output:
        index=expand("out/bwa/index_alt/{{filler}}.{ext}", ext=["amb","ann","bwt","pac","sa"]),
        fasta =      "out/bwa/index_alt/{filler}.fa",
        alt =        "out/bwa/index_alt/{filler}.alt"
    log:
                     "out/bwa/index_alt/{filler}.log"
    benchmark:
                     "out/bwa/index_alt/{filler}.benchmark.tsv"
    params:
        prefix="out/bwa/index_alt/{filler}"
    conda:
        "../envs/bwa.yaml"
    threads:
        1
    shell:
        """
        (ln -srf {input.fasta} {output.fasta}
        ln -srf {input.alt} {output.alt}
        bwa index -p {params.prefix} {output.fasta} ) &> {log}
        """

rule awk_extract_sam_cigar_for_contig_stats:
    """
    Aim:
        Tsv is a subset to make stats on cigar string in R
    Test:
        out/awk/extract_sam_cigar_for_contig_stats/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.tsv
    """
    input:
        sam = "out/{filler}.sam",
    output:
        tsv = "out/awk/extract_sam_cigar_for_contig_stats/{filler}.tsv"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        # Tsv is a subset to make stats on cigar srting in R
        awk 'BEGIN{{FS=OFS="\\t"}}{{if ( $1 ~ !/^@/ ) {{print $3":"$4,$6,$10}} }} ' {input.sam} > {output.tsv}
        """

rule awk_aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3b:
    """
    Aim:
        Take a fasta genome and a sam file that should contain contigs created for example by edena then aligned to the fasta genome/
        Contig with inserts are extracted from the sam file to the output one and their sequence are added to the output fa genome.
    TODO Test with GRCh38 and alt 
    Test:
        out/awk/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_fa-genome-hg19-main-chr/bowtie2/se_fasta_Abraham2017_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_all_DNA_samples_contigs.sam
    """
    input:
        sam="out/{filler}.sam",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        fasta =        "out/awk/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}.fa",
        sam   =        "out/awk/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}.sam",
        alt   =        "out/awk/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}.alt"
    params:
        prefix=        "out/awk/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        # Only keeping insertions strictly greater than 3 to keep alternate contig number low as insertion below this length seems to be correctly called by GATK and Platypus on the reference assembly.
        awk 'BEGIN{{FS=OFS="\\t"}}{{if ($1 ~ /^@/ ) {{print $0}} else if ( $6 ~ /[0-9]{{2,}}I|[4-9]I/ ) {{$1=$3 ":" $4 "_" $6 "_" NR "_alt"; print $0}} }} ' {input.sam} > {output.sam}

        # alt is just the same as sam but this suffix is required by bwa index.
        ln {output.sam} {output.alt}
        
        cat {input.fasta} > {output.fasta}
        awk '{{if ($1 !/^@/ && $6 ~ /I/ ) {{print ">" $3 ":" $4 "_" $6 "_" NR "_alt\\n" $10 }} }} ' {output.sam} >> {output.fasta}
        """

rule awk_aggregate_fa_genome_and_sam_contigs_with_indel_gt_3b:
    """
    Aim:
        Take a fasta genome and a sam file that should contain contigs created for example by edena then aligned to the fasta genome/
        Contig with inserts are extracted from the sam file to the output one and their sequence are added to the output fa genome.
    TODO Test with GRCh38 and alt 
    Test:
        out/awk/aggregate_fa_genome_and_sam_contigs_with_indels_gt_3bp_fa-genome-hg19-main-chr/bowtie2/se_fasta_Abraham2017_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_H3K27ac_contigs.sam
    """
    input:
        sam="out/{filler}.sam",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        fasta =        "out/awk/aggregate_fa_genome_and_sam_contigs_with_indels_gt_3bp_{fa_genome_id}/{filler}.fa",
        sam   =        "out/awk/aggregate_fa_genome_and_sam_contigs_with_indels_gt_3bp_{fa_genome_id}/{filler}.sam",
        alt   =        "out/awk/aggregate_fa_genome_and_sam_contigs_with_indels_gt_3bp_{fa_genome_id}/{filler}.alt"
    params:
        prefix=        "out/awk/aggregate_fa_genome_and_sam_contigs_with_indels_gt_3bp_{fa_genome_id}/{filler}"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        # Only keeping insertions strictly greater than 3 to keep alternate contig number low as insertion below this length seems to be correctly called by GATK and Platypus on the reference assembly.
        awk 'BEGIN{{FS=OFS="\\t"}}{{if ($1 ~ /^@/ ) {{print $0}} else if ( $6 ~ /[0-9]{{2,}}[DI]|[4-9][DI]/ ) {{$1=$3 ":" $4 "_" $6 "_" NR "_alt"; print $0}} }} ' {input.sam} > {output.sam}

        # alt is just the same as sam but this suffix is required by bwa index.
        ln {output.sam} {output.alt}
        
        cat {input.fasta} > {output.fasta}
        awk '{{if ($1 !/^@/ && $6 ~ /[DI]/ ) {{print ">" $3 ":" $4 "_" $6 "_" NR "_alt\\n" $10 }} }} ' {output.sam} >> {output.fasta}
        """



rule meta_aggregate_fa_genome_with_sam_alternative_contigs:
    """
    Aim:
        Take a fasta genome and a sam file that should contain contigs created for example by edena then aligned to the fasta genome/
        Contig with inserts are extracted from the sam file to the output one and their sequence are added to the output fa genome.
        These two files in the output directory can be used to generate a bwa index with alt mapping.
    TODO Test with GRCh38 and alt 
    Test:
        out/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_fa-genome-hg19-main-chr/bowtie2/se_fasta_Abraham2017_hg19/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/T11C_all_DNA_samples_contigs.sam
    """
    input:
        sam="out/{filler}.sam",
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        fasta =        "out/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}.fa",
        index = expand("out/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{{fa_genome_id}}/{{filler}}.{ext}", ext=["amb","ann","bwt","pac","sa"]),
        sam   =        "out/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}.sam",
        alt   =        "out/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}.alt"
    params:
        prefix=        "out/bwa/aggregate_fa_genome_and_sam_contigs_with_inserts_gt_3bp_{fa_genome_id}/{filler}"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        # Only keeping insertions strictly greater than 3 to keep alternate contig number low as insertion below this length seems to be correctly called by GATK and Platypus on the reference assembly.
        awk 'BEGIN{{FS=OFS="\\t"}}{{if ($1 ~ /^@/ ) {{print $0}} else if ( $6 ~ /[^0-9]*[4-9]I/ ) {{$1=$3 ":" $4 "_" $6 "_" NR "_alt"; print $0}} }} ' {input.sam} > {output.sam}

        # alt is just the same as sam but this suffix is required by bwa index.
        ln {output.sam} {output.alt}
        
        cat {input.fasta} > {output.fasta}
        awk '{{if ($1 !/^@/ && $6 ~ /I/ ) {{print ">" $3 ":" $4 "_" $6 "_" NR "_alt\\n" $10 }} }} ' {output.sam} >> {output.fasta}
        
        bwa index -p {params.prefix} {output.fasta}
        """
