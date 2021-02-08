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
        fasta = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
        #fai   = lambda wildcards: [path + '.fai' for path in eval(mwconf['ids'][wildcards.bam_list_id])],
        #genome=input_fa_genome
    output:
        index=expand("out/bwa/index/{{fa_genome_id}}.{ext}", ext=["fa", "amb","ann","bwt","pac","sa"])
    params:
        prefix="out/bwa/index/{fa_genome_id}"
    conda:
        "../envs/bwa.yaml"
    threads:
        1
    shell:
        "ln -srf {input.fasta} out/bwa/index/{wildcards.fa_genome_id}.fa; "
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
        # Tsv is a subset to make stats on cigar string in R
        awk 'BEGIN{{FS=OFS="\\t"}}{{if ( $1 ~ !/^@/ ) {{print $3":"$4,$6,$10}} }} ' {input.sam} > {output.tsv}
        """

rule awk_extract_sam_with_indel:
    """
    Aim:
        Tsv is a subset to make stats on cigar string in R
    Test:
        out/awk/extract_sam_with_indel/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_GRCh38/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac_contigs.fasta
        out/awk/extract_sam_with_indel/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_GRCh38/spades/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/contigs.bed
    """
    input:
        sam = "out/{filler}.sam",
    output:
        fasta = "out/awk/extract_sam_with_indel/{filler}.fasta",
        #fasta_indel_mp10bp = "out/awk/extract_sam_with_indel/{filler}.indel_mp10bp.fasta",
        bed   = "out/awk/extract_sam_with_indel/{filler}.bed"
    conda:
        "../envs/gawk.yaml"
    shell:
        """
        # Tsv is a subset to make stats on cigar string in R
        gawk -v BED={output.bed} -v FASTA={output.fasta} 'BEGIN{{FS=OFS="\\t"}}{{if ( $1 ~ !/^@/ && $6 ~ /[ID]/ ) {{
            CONTIG_CHROMOSOME=$3 ;
            CONTIG_START=$4 ;
            CONTIG_SEQUENCE=$10 ;
            CONTIG_LENGTH=length(CONTIG_SEQUENCE) ;
            CONTIG_END=CONTIG_START+CONTIG_LENGTH ;
            CIGAR=$6 ;
            #CIGAR_NO_MD=CIGAR ;
            #CIGAR_NO_MD=gsub(/[0-9]+[MD]/,"",CIGAR_NO_MD) ;
            #CIGAR_NO_MI=CIGAR ;
            #CIGAR_NO_MI=gsub(/[0-9]+[MI]/,"",CIGAR_NO_MI) ;
            
            # match only get the first insertion occurence,
            # meaning downtream coordinates will be erroneous if more than 1 insertion
            # is in the contig
            # This could help get correct coordinates in such edge cases:
            # https://stackoverflow.com/questions/40569441/awk-match-multiple-matches
            match(CIGAR, /([0-9]+)I/, I_ARRAY)

            # Same for deletion here
            match(CIGAR, /([0-9]+)D/, D_ARRAY)

            CONTIG_END = CONTIG_START + CONTIG_LENGTH + D_ARRAY[1] - I_ARRAY[1]
            
            CIGAR_M_UPSTREAM = CIGAR;
            gsub(/M.*$/,"",CIGAR_M_UPSTREAM)
            CIGAR_M_DOWNSTREAM = CIGAR;
            gsub(/^.*[ID]/,"",CIGAR_M_DOWNSTREAM)
            gsub(/M/,"",CIGAR_M_DOWNSTREAM)

            if (0+CIGAR_M_DOWNSTREAM < 10) {{
                INDEL_PLUS_10BP_END=CONTIG_END
                INDEL_MP_10BP_SEQUENCE=CONTIG_SEQUENCE
            }} else {{
                INDEL_PLUS_10BP_END = CONTIG_END + 10 - CIGAR_M_DOWNSTREAM;
                INDEL_MP_10BP_SEQUENCE=substr(CONTIG_SEQUENCE, 1, length(CONTIG_SEQUENCE) - CIGAR_M_DOWNSTREAM + 10)
                }}

            if (0+CIGAR_M_UPSTREAM < 10) {{
                INDEL_MINUS_10BP_START = CONTIG_START
            }} else {{
                INDEL_MINUS_10BP_START = CONTIG_START - 10 + CIGAR_M_UPSTREAM
                INDEL_MP_10BP_SEQUENCE=substr(INDEL_MP_10BP_SEQUENCE, CIGAR_M_UPSTREAM - 10)
                }}

            N_MATCHED_UPSTREAM=CIGAR ;
            sub(/M.*$/,"",N_MATCHED_UPSTREAM) ;
            N_MATCHED_DOWNSTREAM=CIGAR ;
            sub(/M$/,"",N_MATCHED_DOWNSTREAM) ;
            sub(/^.*[ID]/,"",N_MATCHED_DOWNSTREAM) ;
            #INDEL_SEQUENCE=
            #INDEL_NAME=CONTIG_CHROMOSOME":"CONTIG_START"-"CONTIG_END"_"CIGAR"_"NR
            INDEL_NAME=CONTIG_CHROMOSOME":"INDEL_MINUS_10BP_START"-"INDEL_PLUS_10BP_END"_"CIGAR"_"NR

            #print CONTIG_CHROMOSOME,CONTIG_START,CONTIG_END,INDEL_NAME,"0","." > BED
            
            #print CONTIG_CHROMOSOME,CONTIG_START,CONTIG_END > BED
            print CONTIG_CHROMOSOME,INDEL_MINUS_10BP_START,INDEL_PLUS_10BP_END > BED

            print ">" INDEL_NAME > FASTA
            print INDEL_MP_10BP_SEQUENCE > FASTA
            #print CONTIG_SEQUENCE > FASTA
            
            }} }}' {input.sam}
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
        fasta = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
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
        fasta = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
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
        fasta = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
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
