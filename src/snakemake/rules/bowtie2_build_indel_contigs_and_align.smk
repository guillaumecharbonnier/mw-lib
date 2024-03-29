rule bowtie2_build_indel_contigs_and_align_se:
    """
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
        out/awk/extract_sam_with_indel/bowtie2/se_fasta_Abraham2017_GRCh38/spades/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/contigs.bed
        out/bowtie2_build_indel_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.sam
    """
    input:
        #fasta="out/{filler}",
        fastq = "out/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.fastq.gz",
        fasta_indel = "out/awk/extract_sam_with_indel/bowtie2/se_fasta_Abraham2017_GRCh38/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac_contigs.fasta",
        fasta_reference = "out/bedtools/getfasta_fa-genome-hg38/sort/_-u/awk/extract_sam_with_indel/bowtie2/se_fasta_Abraham2017_GRCh38/edena/assembling_-d_20_-c_20_-minCoverage_5/edena/overlapping/gunzip/to-stdout/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac_contigs.fa"
    output:
        sam = "out/bowtie2_build_indel_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.sam",
        unmapped = "out/bowtie2_build_indel_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.fastq.gz",
        bt1 = "out/bowtie2_build_indel_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.1.bt2"
    log:
        "out/bowtie2_build_indel_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.log"
    params:
        index="out/bowtie2_build_indel_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac"
        #index="out/bowtie2-build/{filler}"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        MAX_THREADS
    shell:
        """
        ln -srf {input.fasta_reference} {params.index}.reference.fa
        ln -srf {input.fasta_indel} {params.index}.indel.fa
        bowtie2-build {params.index}.reference.fa,{params.index}.indel.fa {params.index}

        bowtie2 -p {threads} -x {params.index}\
            --very-sensitive-local\
            -U {input.fastq} \
            --un-gz {output.unmapped} \
            -S {output.sam}    2> {log}
        """

rule bowtie2_build_indel_spades_contigs_and_align_se:
    """
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Test:
        out/samtools/idxstats/samtools/index/picard/MarkDuplicates_REMOVE_DUPLICATES=true_VALIDATION_STRINGENCY=SILENT/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2_build_indel_spades_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.idxstat.tsv
    """
    input:
        #fasta="out/{filler}",
        fastq = "out/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.fastq.gz",
        fasta_indel = "out/awk/extract_sam_with_indel/bowtie2/se_fasta_Abraham2017_GRCh38-ensembl-r100/spades/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/contigs.fasta",
        fasta_reference = "out/bedtools/getfasta_fa-genome-GRCh38-ensembl-r100/sort/_-u/awk/extract_sam_with_indel/bowtie2/se_fasta_Abraham2017_GRCh38-ensembl-r100/spades/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/contigs.fa"
    output:
        sam      = "out/bowtie2_build_indel_spades_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.sam",
        unmapped = "out/bowtie2_build_indel_spades_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.fastq.gz",
        bt1      = "out/bowtie2_build_indel_spades_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.1.bt2"
    log:
                   "out/bowtie2_build_indel_spades_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.log"
    params:
        index    = "out/bowtie2_build_indel_spades_contigs_and_align_se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac"
        #index="out/bowtie2-build/{filler}"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        MAX_THREADS
    shell:
        """
        ln -srf {input.fasta_reference} {params.index}.reference.fa
        ln -srf {input.fasta_indel} {params.index}.indel.fa
        bowtie2-build {params.index}.reference.fa,{params.index}.indel.fa {params.index}

        bowtie2 -p {threads} -x {params.index}\
            --very-sensitive-local\
            -U {input.fastq} \
            --un-gz {output.unmapped} \
            -S {output.sam} #| samtools view -L {input.bed_indel_coords} 2> {log}
        """

rule bowtie2_build_flanked_indel_spades_and_align_se:
    """
    Created:
        2020-08-22 15:31:58
    Test:
        out/bowtie2/build_flanked_indel_spades_and_align_se_GRCh38-ensembl-r100/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.sam
    """
    input:
        fasta="out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/grep/extract-indel-from-vcf/bcftools/mpileup_fa-genome-{genome_id}/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_Abraham2017_{genome_id}/spades/se/{filler}/contigs.fasta",
        fastq="out/{filler}.fastq.gz"
    output:
        sam      = "out/bowtie2/build_flanked_indel_spades_and_align_se_{genome_id}/{filler}.sam",
        unmapped = "out/bowtie2_build_flanked_indel_spades_and_align_se_{genome_id}/{filler}.fastq.gz",
    log:
        "out/bowtie2/build_flanked_indel_spades_and_align_se_{genome_id}/{filler}.log"
    params:
        index = "out/bowtie2/build_flanked_indel_spades_and_align_se_{genome_id}/{filler}"
    wildcard_constraints:
        genome_id = "[a-zA-Z0-9-]+"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        (
        # Softlink here because bowtie2-build
        # misunderstands the comma in input.fasta
        # for a separator
        ln -srf {input.fasta} {params.index}.fa

        bowtie2-build {params.index}.fa {params.index}

        bowtie2 -p {threads} -x {params.index}\
            --very-sensitive-local\
            -U {input.fastq} \
            --un-gz {output.unmapped} \
            -S {output.sam}
        
        ) 2> {log}
        """


rule bowtie2_build_flanked_indel_bbnorm_spades_and_align_se:
    """
    Created:
        2020-08-22 15:31:58
    Test:
        out/bowtie2/build_flanked_indel_bbnorm_spades_and_align_se_GRCh38-ensembl-r100/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.sam
    """
    input:
        fasta="out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/grep/extract-indel-from-vcf/bcftools/mpileup_fa-genome-{genome_id}/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_Abraham2017_{genome_id}/spades/se/bbnorm/se/{filler}/contigs.fasta",
        fastq="out/{filler}.fastq.gz"
    output:
        sam      = "out/bowtie2/build_flanked_indel_bbnorm_spades_and_align_se_{genome_id}/{filler}.sam",
        unmapped = "out/bowtie2_build_flanked_indel_bbnorm_spades_and_align_se_{genome_id}/{filler}.fastq.gz",
    log:
        "out/bowtie2/build_flanked_indel_bbnorm_spades_and_align_se_{genome_id}/{filler}.log"
    params:
        index = "out/bowtie2/build_flanked_indel_bbnorm_spades_and_align_se_{genome_id}/{filler}"
    wildcard_constraints:
        genome_id = "[a-zA-Z0-9-]+"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        (
        # Softlink here because bowtie2-build
        # misunderstands the comma in input.fasta
        # for a separator
        ln -srf {input.fasta} {params.index}.fa

        bowtie2-build {params.index}.fa {params.index}

        bowtie2 -p {threads} -x {params.index}\
            --very-sensitive-local\
            -U {input.fastq} \
            --un-gz {output.unmapped} \
            -S {output.sam}
        
        ) 2> {log}
        """

rule bowtie2_build_flanked_indel_megahit_and_align_se:
    """
    Created:
        2020-08-22 15:31:58
    Test:
        out/bowtie2/build_flanked_indel_megahit_and_align_se_GRCh38-ensembl-r100/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac.sam
    """
    input:
        # Before 2023-07-06, this was:
        # fasta="out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/python/extract_indels_from_sam_to_vcf_fa-genome-{genome_id}/bowtie2/se_fa_Abraham2017_{genome_id}/megahit/se{megahit_args}/{filler}.contigs.fasta",
        # After 2023-07-06, we only keep variants that are not in dbSNP (i.e. novel variants that could cause neoenhancer):
        fasta = "out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/bcftools/view_vcf_to_vcf_keep-no-id/snpsift/annotate_dbsnp/python/extract_indels_from_sam_to_vcf_fa-genome-GRCh38-ensembl-r100/bowtie2/se_fa_Abraham2017_GRCh38-ensembl-r100/megahit/se{megahit_args}/{filler}.contigs.fasta",
        fastq="out/{filler}.fastq.gz"
    output:
        sam      = "out/bowtie2/build_flanked_indel_megahit{megahit_args}_and_align_se_{genome_id}/{filler}.sam"
    log:
        "out/bowtie2/build_flanked_indel_megahit{megahit_args}_and_align_se_{genome_id}/{filler}.log"
    params:
        index = "out/bowtie2/build_flanked_indel_megahit{megahit_args}_and_align_se_{genome_id}/{filler}"
    wildcard_constraints:
        genome_id = "[a-zA-Z0-9-]+",
        megahit_args = "[a-zA-Z0-9-_]*"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        (
        # Softlink here because bowtie2-build
        # misunderstands the comma in input.fasta
        # for a separator
        ln -srf {input.fasta} {params.index}.fa

        bowtie2-build {params.index}.fa {params.index}

        bowtie2 -p {threads} -x {params.index}\
            --very-sensitive-local\
            -U {input.fastq} \
            --no-unal -S {output.sam}

        # We can the remove the index because it is specific of this sample
        rm {params.index}*.bt2
        
        ) 2> {log}
        """

rule bowtie2_build_flanked_indel_megahit_and_align_pe:
    """
    Created:
        2020-08-23 18:00:26
    Note:
        I add {megahit_args} on 2023-06-17 because I am not sure about the most appropriate settings between default and --min-count 3 --no-mercy.
        The uncertainty arose because I find only a few variants using the --min-count 3 --no-mercy settings.
        The code is adapted to be able to generate html report for any megahit_args in order to investigate this.
    Test:
        out/bowtie2/build_flanked_indel_megahit_and_align_pe_GRCh38-ensembl-r100/ln/alias/sst/all_samples/fastq/856_H3K27ac.sam

        out/python/vcf_to_flanking_sequence_fa-genome-GRCh38-ensembl-r100/grep/extract-indel-from-vcf/bcftools/mpileup_fa-genome-GRCh38-ensembl-r100/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fa_Abraham2017_GRCh38-ensembl-r100/megahit/pe/ ln/alias/sst/all_samples/fastq/856_H3K27ac/final.contigs.fasta
    """
    input:
        # # fasta="out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/grep/extract-indel-from-vcf/bcftools/mpileup_fa-genome-{genome_id}/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fa_Abraham2017_{genome_id}/megahit/pe_--min-count_3_--no-mercy/{filler}/final.contigs.fasta",
        # # Trying to revert to the default megahit settings because there is not many calls with no-mercy
        # fasta="out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/python/extract_indels_from_sam_to_vcf_fa-genome-{genome_id}/bowtie2/se_fa_Abraham2017_{genome_id}/megahit/pe{megahit_args}/{filler}.contigs.fasta",
        # Before 2023-07-06, this was:
        # fasta="out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/python/extract_indels_from_sam_to_vcf_fa-genome-{genome_id}/bowtie2/se_fa_Abraham2017_{genome_id}/megahit/pe{megahit_args}/{filler}.contigs.fasta",
        # After 2023-07-06, we only keep variants that are not in dbSNP (i.e. novel variants that could cause neoenhancer):
        fasta = "out/python/vcf_to_flanking_sequence_fa-genome-{genome_id}/bcftools/view_vcf_to_vcf_keep-no-id/snpsift/annotate_dbsnp/python/extract_indels_from_sam_to_vcf_fa-genome-GRCh38-ensembl-r100/bowtie2/se_fa_Abraham2017_GRCh38-ensembl-r100/megahit/pe{megahit_args}/{filler}.contigs.fasta",
        fastq_1="out/{filler}_1.fastq.gz",
        fastq_2="out/{filler}_2.fastq.gz"
    output:
        sam      = "out/bowtie2/build_flanked_indel_megahit{megahit_args}_and_align_pe_{genome_id}/{filler}.sam"
    log:
        "out/bowtie2/build_flanked_indel_megahit{megahit_args}_and_align_pe_{genome_id}/{filler}.log"
    params:
        index = "out/bowtie2/build_flanked_indel_megahit{megahit_args}_and_align_pe_{genome_id}/{filler}"
    wildcard_constraints:
        genome_id = "[a-zA-Z0-9-]+",
        megahit_args = "[a-zA-Z0-9-_]*"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        (
        # Softlink here because bowtie2-build
        # misunderstands the comma in input.fasta
        # for a separator
        ln -srf {input.fasta} {params.index}.fa

        bowtie2-build {params.index}.fa {params.index}

        bowtie2 -p {threads} -x {params.index}\
            --very-sensitive-local\
            -1 {input.fastq_1} -2 {input.fastq_2}\
            --no-unal -S {output.sam}

        # We can the remove the index because it is specific of this sample
        rm {params.index}*.bt2
        
        ) 2> {log}
        """



