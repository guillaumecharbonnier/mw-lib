rule sed_extra:
    """
    Created:
        2019-03-15 08:43:31
    Test:
        out/sed/rm-chr-in-fa/cat/assembly_ensembl/NCBIM37.fa
        Fastq from SOLID and converted by TGML workflow have a ! at the begining of qual line.
        I can remove it with this rule to match fastq produced by xsqtools+bbfast solid2fastq:
        out/sed/rm-first-exclam-in-fq/ln/abspath/gpfs/tagc/home/sadouni/reads/fastq/Run_75/Library_1/Run_75_L01_Library_1_F3.fastq
    """
    input:
        "out/{filler}"
    output:
        temp("out/{tool}{extra}/{filler}")
    params:
        extra = params_extra
    conda:
        "../envs/coreutils.yaml"
    wildcard_constraints:
        tool="sed/"
    shell:
        "sed {params} {input} > {output}"


ruleorder: sed_add_chr_legacy > sed_extra
rule sed_add_chr_legacy:
    """
    Created:
        2017-03-08 14:18:47
    Aim:
        gtf from Ensembl use 1,2,3... instead of chr1,chr2,chr3... and this is annoying so we add chr before to have the same chr names that in chrominfo and elsewhere.
    Test:
        out/sed/add_chr/gunzip/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
    """
    input:
        "out/{filler}"
    output:
        "out/sed/add_chr/{filler}"
    shell:
        """
        sed '/^#/! s/^/chr/g' {input} > {output}
        """

ruleorder: sed_add_chr_fa_legacy > sed_extra
rule sed_add_chr_fa_legacy:
    """
    Created:
        2018-02-21 15:09:10
    Aim:
        Add 'chr' to chromosome name from Ensembl to make them like in UCSC.
    Test:
        out/sed/add_chr_fa/gunzip/to-stdout/wget/ftp_ensembl/pub/release-67/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.67.dna.chromosome.1.fa
    """
    input:
        "out/{filler}"
    output:
        "out/sed/add_chr_fa/{filler}"
    shell:
        """
        sed '/^#/! s/^>/>chr/g' {input} > {output}
        """

ruleorder: sed_rename_mitochondrial_chr_M_to_MT_legacy > sed_extra
rule sed_rename_mitochondrial_chr_M_to_MT_legacy:
    """
    Created:
        2017-05-16 16:53:58
    Aim:
        Because in mm10 the mitochondrial chromosome is 'M' whereas it is 'MT' in GRCm38...
    Test:
        out/sed/rename_mitochondrial_chr_M_to_MT/sed/remove_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt
    """
    input:
        "out/{filler}"
    output:
        "out/sed/rename_mitochondrial_chr_M_to_MT/{filler}"
    shell:
        """
        sed 's/^M\t/MT\t/g' {input} > {output}
        """

ruleorder: sed_remove_chr_legacy > sed_extra
rule sed_remove_chr_legacy:
    """
    Created:
        2017-03-28 16:09:54
    Aim:
        When working with ENSEMBL reference the genomes do not have the 'chr' prefix and this can cause issue when using a reference having them like Rseqc housekeeping bed.
    Test:
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/sed/remove_chr/{filler}"
    shell:
        """
        sed 's/^chr//' {input.txt} > {output.txt}
        """

ruleorder: sed_remove_m_rand_un_chr_legacy > sed_extra
rule sed_remove_m_rand_un_chr_legacy:
    """
    Created:
        2017-11-12 23:17:03
    Aim:
    Test:
        out/sed/remove_m_rand_un_chr/bowtie2/pe_mm9/sickle/pe_-t_sanger_-q_30/ln/paired_end_remove_mate_prefix/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run140/Input-Nut-WT.sam
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/sed/remove_m_rand_un_chr/{filler}"
    shell:
        """
        sed '/chrM/d;/random/d;/chrUn/d' {input.txt} > {output.txt}
        """

ruleorder: sed_mac_to_unix_legacy > sed_extra
rule sed_mac_to_unix_legacy:
    """
    Created:
        2017-10-02 14:02:11
    Aim:
        Some files from Sophie Rousseaux come with weird end of line (^M character displayed but not taken correctly by unix). Replacing it with newlines does the trick.
    Note:
        The command is inside a script because snakemake shell does'not like the ^M.
    Test:
       out/sed/mac_to_unix/inp/sophie_rousseaux_2017_09_29/anova_nut_sp_line_res_3_RSkovRSwt_gsea1.txt 
    """
    input:
        txt="out/{filler}",
        src="src/bash/sed_mac_to_unix.sh"
    output:
        txt="out/sed/mac_to_unix/{filler}"
    shell:
        """
        {input.src} {input.txt} > {output.txt}
        """

ruleorder: sed_gene_name_to_uppercase_in_gtf_legacy > sed_extra
rule sed_gene_name_to_uppercase_in_gtf_legacy:
    """
    Created:
        2017-10-02 15:32:55
    Aim:
        Convert gene names in gtf from mixed case to uppercase. This is done in order to be able to retrieve gene_name from nut microarray which are all uppercase.
    Test:
        input:
            out/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
        output:
            out/sed/gene_name_to_uppercase_in_gtf/gunzip/to-stdout/wget/ftp_ensembl/pub/release-87/gtf/mus_musculus/Mus_musculus.GRCm38.87.gtf
    """
    input:
        gtf="out/{filler}.gtf"
    output:
        gtf="out/sed/gene_name_to_uppercase_in_gtf/{filler}.gtf"
    shell:
        """
        sed 's/gene_name "\([a-zA-Z0-9]*\)"/gene_name "\\U\\1"/g' {input.gtf} > {output.gtf}
        """

ruleorder: sed_remove_first_and_last_line_legacy > sed_extra
rule sed_remove_first_and_last_line_legacy:
    """
    Created:
        2018-03-13 16:47:29
    Aim:
        Because GREAT Txt files add a header and an empty last line.
    Test:
        input:
            out/doc/great/hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.txt
        output:
            out/sed/remove_first_and_last_line/doc/great/hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.txt
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/sed/remove_first_and_last_line/{filler}"
    shell:
        """
        sed '1d;$d' {input.txt} > {output.txt}
        """

ruleorder: sed_extract_go_id_from_obo_legacy > sed_extra
rule sed_extract_go_id_from_obo_legacy:
    """
    Created:
        2018-11-23 20:38:27
    Test:
        out/sed/extract_go_id_from_obo/wget/http/www.geneontology.org/ontology/subsets/goslim_generic.txt
    """
    input:
        "out/{filler}.obo"
    output:
        "out/sed/extract_go_id_from_obo/{filler}.txt"
    shell:
        "sed -n 's/^.*id: GO/GO/p' {input} > {output}"

ruleorder: sed_extract_great_go_biological_process_legacy > sed_extra
rule sed_extract_great_go_biological_process_legacy:
    """
    Created:
        2018-04-12 16:17:52
    Aim:
        GREAT all.tsv are painful to read directly in R because there are some '#' inside MsigDB terms. I do not need them so I just extract first my rows of interest before importing to R.
    """
    input:
        tsv="out/{filler}"
    output:
        tsv="out/sed/extract_great_go_biological_process/{filler}"
    shell:
        "sed -n '/^GO Biological Process/p' {input.tsv} > {output.tsv}"

rule r_heatmap_from_great_top_go_biological_process:
    """
    Created:
        2018-05-15 22:58:41
    Test:
        out/r/heatmap_from_great_top-10_go_biological_process/heatmap.pdf out/r/heatmap_from_great_top-20_go_biological_process/heatmap.pdf
    """
    input:
        Rscript="opt/miniconda/envs/r/bin/Rscript",
        script="src/r/script/heatmap_from_great_top_go_biological_process.R",
        great_process=expand("out/sed/extract_great_go_biological_process/doc/great/hg19/crossmap/bed_hg38_to_hg19/awk/fill_bed3_to_bed6/bedtools/multiinter_thymus_peaks_hg38_merged/{sample}.tsv",sample=["1","2","3","4","5"])
    output:
        pdf="out/r/heatmap_from_great_top-{top}_go_biological_process/heatmap.pdf"
    shell:
        "{input.Rscript} {input.script} -p {output.pdf} -t {wildcards.top} {input.great_process}"


rule r_heatmap_from_great_top_20_go_biological_process:
    """
    Created:
        2018-05-15 22:58:41
    """
    input:
        Rscript="opt/miniconda/envs/r/bin/Rscript",
        script="src/r/script/heatmap_from_great_top_go_biological_process.R",
        great_process=expand("out/sed/extract_great_go_biological_process/doc/great/hg19/crossmap/bed_hg38_to_hg19/awk/fill_bed3_to_bed6/bedtools/multiinter_thymus_peaks_hg38_merged/{sample}.tsv",sample=["1","2","3","4","5"])
    output:
        pdf="out/r/heatmap_from_great_top_go_biological_process/heatmap.pdf"
    shell:
        "{input.Rscript} {input.script} -p {output.pdf} {input.great_process}"

rule r_heatmap_from_great_top_go_biological_process_proximal_and_distal:
    input:
        Rscript="opt/miniconda/envs/r/bin/Rscript",
        script="src/r/script/heatmap_from_great_top_go_biological_process.R",
        great_process=expand("out/sed/extract_great_go_biological_process/doc/great/hg19/awk/extract_distal_proximal/paste/great_input_and_output_hg19/crossmap/bed_hg38_to_hg19/awk/fill_bed3_to_bed6/bedtools/multiinter_thymus_peaks_hg38_merged/{sample}_{distal_proximal}_all.tsv",sample=["1","2","2,3","3","4","4,5","5"], distal_proximal=["proximal","distal"])
    output:
        pdf="out/r/heatmap_from_great_top_go_biological_process_proximal_distal/heatmap.pdf"
    shell:
        "{input.Rscript} {input.script} -p {output.pdf} {input.great_process}"

rule r_heatmap_from_great_top_go_biological_process_proximal_or_distal:
    input:
        Rscript="opt/miniconda/envs/r/bin/Rscript",
        script="src/r/script/heatmap_from_great_top_go_biological_process.R",
        great_process=expand("out/sed/extract_great_go_biological_process/doc/great/hg19/awk/extract_distal_proximal/paste/great_input_and_output_hg19/crossmap/bed_hg38_to_hg19/awk/fill_bed3_to_bed6/bedtools/multiinter_thymus_peaks_hg38_merged/{sample}_{distal_proximal}_all.tsv",sample=["1","2","2,3","3","4","4,5","5"], distal_proximal=["distal"])
    output:
        pdf="out/r/heatmap_from_great_top_go_biological_process_distal/heatmap.pdf"
    shell:
        "{input.Rscript} {input.script} -p {output.pdf} {input.great_process}"

