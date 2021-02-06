rule gfold_diff:
    """
    Created:
        2017-03-15 16:31:53
    Aim:
        Gfold rule using tsv as suffix instead of the uselessly explicit read_cnt.
    Test:
        out/gfold/diff/gfold/count_bed-h4k5ac-peaks/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.tsv
        out/gfold/diff/gfold/count_bed-mm10-rmsk/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-WT_VS_H4K5ac-Nut-KO.tsv
        
        out/gfold/diff/gfold/count_gtf-GRCm38/star/pe_GRCm38/sickle/pe_-t_sanger_-q_20/gunzip/nut_rnaseq/Nut-R-WT_VS_Nut-R-KO.tsv
        TODO: MAKE NUT UPREGDOWNREG WORK AGAIN WITH MM10
    """
    input:
        tsv1="out/{filler}/{sample1}.tsv",
        tsv2="out/{filler}/{sample2}.tsv"
    output:
        tsv="out/gfold/diff/{filler}/{sample1}_VS_{sample2}.tsv"
    params:
        unsuffixed_tsv1="out/{filler}/{sample1}",
        unsuffixed_tsv2="out/{filler}/{sample2}"
    wildcard_constraints:
        sample1="[a-zA-Z0-9_-]+",
        sample2="[a-zA-Z0-9_-]+"
    conda:
        "../envs/gfold.yaml"
    shell:
        """
        gfold diff \
            -s1 {params.unsuffixed_tsv1} \
            -s2 {params.unsuffixed_tsv2} \
            -suf .tsv -o {output.tsv}
        """

rule gfold_diff_with_replicates:
    """
    Created:
        2017-03-27 10:09:38
    Aim:
        Gfold taking replicate list as instead of just one sample per condition.
    Test:
        out/gfold/diff_with_replicates/tsv-gfold-GRCh38-casero2016-thy3_VS_tsv-gfold-GRCh38-casero2016-thy4.tsv
    """
    input:
        tsv1 = lambda wildcards: eval(mwconf['ids'][wildcards.tsv1_id]),
        tsv2 = lambda wildcards: eval(mwconf['ids'][wildcards.tsv2_id]),
    output:
        tsv="out/gfold/diff_with_replicates/{tsv1_id}_VS_{tsv2_id}.tsv"
    wildcard_constraints:
        tsv1_id="tsv-gfold-[a-zA-Z0-9_-]+",
        tsv2_id="tsv-gfold-[a-zA-Z0-9_-]+"
    conda:
        "../envs/gfold.yaml"
    shell:
        """
        # Remove suffix and replace spaces by commas.
        tsv1_list=`echo '{input.tsv1}' | sed 's/\.tsv//g' | tr ' ' ','`
        tsv2_list=`echo '{input.tsv2}' | sed 's/\.tsv//g' | tr ' ' ','`

        echo $tsv1_list

        gfold diff \
            -s1 $tsv1_list \
            -s2 $tsv2_list \
            -suf .tsv \
            -o {output.tsv}
        """

