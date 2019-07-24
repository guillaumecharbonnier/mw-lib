rule coreutils_format_list_of_rules:
    output:
        "out/coreutils/format_list_of_rules/mw.txt"
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        ls -1 ../mw/src/snakemake/rules/ | sed 's/.smk//' > {output}
        """

rule coreutiles_yaml_to_gprofiler_list:
    """
    Created:
        2018-11-22 00:12:45
    Test:
        out/r/chen_sup_to_broad_or_sharp_tables/broad.yaml
        out/coreutils/yaml_to_gprofiler_list/r/chen_sup_to_broad_or_sharp_tables/broad.txt

        cat out/r/chen_sup_to_broad_or_sharp_tables/broad.yaml | sed  's/^[^-]/> \1/'

        http://metascape.org/gp/index.html#/main/step1
    """
    input:
        yaml="out/{filler}.yaml"
    output:
        txt="out/coreutils/yaml_to_gprofiler_list/{filler}.txt",
        tsv="out/coreutils/yaml_to_gprofiler_list/{filler}.tsv"
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        cat {input.yaml} | sed  's/^\([^-]\)/> \\1/' | tr '\\n' ' ' | sed -e 's/>/\\n>/g' -e 's/:/:\\n/g' -e 's/ - / /g' > {output.txt}
        cat {input.yaml} | sed 's/^\([^-]\)/> \\1/' | tr '\\n' ' ' | sed -e 's/> /\\n/g' -e 's/:/\\t/g' -e 's/ - / /g' > {output.tsv}
        """

rule coreutils_extract_min_dist_from_great:
    """
    Created:
        2018-03-13 18:13:27
    Aim:
        GREAT produces output like this:
        unnamed SAMD11 (-147061), OR4F16 (-92004)
        unnamed OR4F16 (-140629)
        unnamed OR4F16 (-151691), SAMD11 (-87374)
        and I want to extract the min value into a one column file to then paste with the original input bed file.
    Test:
        input:
            out/sed/remove_first_and_last_line/doc/great/hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.txt
        output:
            out/coreutils/extract_min_dist_from_great/sed/remove_first_and_last_line/doc/great/hg19/crossmap/bed_hg38_to_hg19/bedtools/intersect_thymus_stage_common_peaks_hg38_v3/not_common_subset_3.txt
    """
    input:
        txt="out/{filler}.txt"
    output:
        txt="out/coreutils/extract_min_dist_from_great/{filler}.txt"
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        rm -f {output.txt}

        while read l;
        do
            echo $l |\
            sed -e 's/(-/;/g' -e 's/)/;/g' -e 's/(+/;/g' |\
            cut -f2,4 -d ';' |\
            tr ';' '\n' |\
            sort -n |\
            head -1 >> {output.txt}
        done < {input.txt}
        """

rule coreutils_join_cpg_matrix:
    input:
        txt=expand("out/gunzip/to-stdout/ln/alias/experiments/thymus_BS/{sample}.txt", sample=BLUEPRINT_THYMUS_BS_CPG)
    output:
        txt="out/coreutils/join_cpg_matrix.txt"
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        FIRST_TXT=`echo '{input.txt}' | cut -f1 -d ' '`

        head $FIRST_TXT | cut -f1,2,6 | sed 's/\\t/-/' > {output.txt}

        for TXT in `echo {input.txt} | cut -f2- -d ' '`;
        do
            cut -f1,2,6 $TXT | sed 's/\\t/-/' | head | join {output.txt} -
        done
        """

rule coreutils_subsample:
    """
    Created: 2016-12-20 15h21 - Created to do random sampling in a bed file but could be applied to any text file.

    Note:
    Superslow compared to shuf, prefer shuf_subsample rule.

    wc -l out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.bed
    69771672 out/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-KO.bed

    The number of line in this file should be used to subsample PSK-SC-WT

    Test:
        "out/sort/subsample69771672/awk/filter_bedpe_by_size/min40_max80/samtools/view/bampe_to_bed3_nodup/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT.bed"
    """
    input:
        txt="out/{filler}"
    output:
        txt="out/coreutils/subsample{n}/{filler}"
        #txt="out/sort/subsample{n}/{filler}"
    wildcard_constraints:
        n="[0-9]+"
    shell:
        """
        sort --random-sort {input.txt} | head -n {wildcards.n} > {output.txt}
        """
