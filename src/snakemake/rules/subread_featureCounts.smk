rule subread_featureCounts_extra:
    """
    Created:
        2018-10-24 20:26:29
    Modified:
        2019-03-03 20:57:04 - Generalization to work with both gtf and saf.
    Aim:
        Counts
    Doc:
        http://gensoft.pasteur.fr/docs/subread/1.4.6-p3/SubreadUsersGuide.pdf#page=27
    Notes from doc:
        -t <string>     Specify feature type in GTF annotation. `exon' by 
                        default. Features used for read counting will be 
                        extracted from annotation using the provided value.

        -g <string>     Specify attribute type in GTF annotation. `gene_id' by 
                        default. Meta-features used for read counting will be 
                        extracted from annotation using the provided value.
    Warning:
        -O argument has to be used for transcript quantification as if a gene has multiple overlapping transcripts, reads shared by multiple genes will not be counted without. Maybe also for gene quantification...
    Tests:
        out/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-GRCh38-ensembl_bam-hg38-Casero2016-thy3-thy4.tsv
        out/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-GRCh38-ensembl_bam-hg38-RNA-thymus.tsv
        out/subread/featureCounts_-F_SAF_saf-mm10-test-macs2-peaks_bam-mm10-test-srr.tsv
        out/subread/featureCounts_-F_SAF_saf-mm10-tss-pm1kb_bam-mnase-cs-spm-nuc-ss.tsv
    """
    input:
        gtf_saf = lambda wildcards: config['ids'][wildcards.gtf_saf_id],
        bam     = lambda wildcards: eval(config['ids'][wildcards.bam_list_id])
    output:
        tsv = "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.tsv",
        txt = "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.tsv.summary"
    log:
              "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.log"
    benchmark:
              "out/{tool}{extra}_{gtf_saf_id}_{bam_list_id}.benchmark.tsv"
    params:
        extra = params_extra
    conda:
        "../envs/subread.yaml"
    wildcard_constraints:
        tool="subread/featureCounts",
        bam_list="bam-[a-zA-Z0-9-]+",
        gtf_saf_id="[a-zA-Z0-9-]+"
    threads:
        16
    shell:
        """
        featureCounts -a {input.gtf_saf} {params.extra} -T {threads} -o {output.tsv} {input.bam}
        """



#
#rule featureCounts:
#    """
#    Created:
#        2016-11-21 10h17 
#    Notes from doc:
#        -t <string>     Specify feature type in GTF annotation. `exon' by 
#                        default. Features used for read counting will be 
#                        extracted from annotation using the provided value.
#
#        -g <string>     Specify attribute type in GTF annotation. `gene_id' by 
#                        default. Meta-features used for read counting will be 
#                        extracted from annotation using the provided value.
#
#    Obsolete tests:
#        "out/featureCounts/gtfGRCm38_texon_ggene_id_test.tsv"
#        "out/featureCounts/gtfGRCm38_tfive_prime_utr_ggene_id_test.tsv"
#        "out/featureCounts/gtfGRCm38_texon_ggene_id_khan_2017_02_14.tsv"
#        out/featureCounts/gtfmm9CustomFt_ttest_ggene_id_nut_chipseq.tsv
#        out/featureCounts/gtfmerge-attr-GRCm38_texon_gmerge_gene_id_name_h2al2_rnaseq.tsv
#        out/featureCounts/gtfcuffmerge-GRCm38-H2AL2_texon_ggene_id_h2al2_rnaseq.tsv
#        out/featureCounts/gtfcuffmerge-GRCm38-H2AL2-outFilterMultimapNmax-1000_texon_ggene_id_h2al2_rnaseq.tsv
#        out/featureCounts/gtf-longs-merge-attr-GRCm38_t-exon_g-merge_gene_id_name_bam-h2al2-rnaseq.tsv
#        out/featureCounts/gtf-longs-protein-coding-merge-attr-GRCm38_t-exon_g-merge_gene_id_name_bam-h2al2-rnaseq.tsv
#    Tests:
#        out/subread/featureCounts_gtf-hg38-ensembl_t-exon_g-transcript_id_bam-hg38-RNA-thymus.tsv
#        out/subread/featureCounts_O_gtf-hg38-ensembl_t-exon_g-transcript_id_bam-hg38-RNA-thymus.tsv
#        out/subread/featureCounts_M_p_gtf-hg38-ensembl_t-exon_g-transcript_id_bam-hg38-RNA-thymus.tsv
#    Warning:
#        -O argument has to be used for transcript quantification as if a gene has multiple overlapping transcripts, reads shared by multiple genes will not be counted without. Maybe also for gene quantification...
#        -ignoreDup
#    """
#    input:
#        gtf=input_gtf,
#        bam=input_bam_list
#    output:
#        tsv="out/subread/featureCounts{extra}_gtf-{gtf_id}_t-{feature_type}_g-{attribute_type}_bam-{bam_list_id}.tsv"
#    conda:
#        "../envs/subread.yaml"
#    wildcard_constraints:
#        feature_type="exon|five_prime_utr|test",
#        attribute_type="gene_id|merge_gene_id_name|transcript_id",
#        bam_list="[a-zA-Z0-9-]+"
#    threads:
#        1
#        #16 Even with more threads featureCounts only use one...
#    shell:
#        """
#        EXTRA=`echo {wildcards.extra} | sed -e 's/-/ /g' -e 's/_/ -/g'`
#        featureCounts -a {input.gtf} -t {wildcards.feature_type} -g {wildcards.attribute_type} $EXTRA -T {threads} -o {output.tsv} {input.bam}
#        """
#
rule awk_extract_counts_from_featureCounts:
    """
    Created:
        2017-02-15 11:52:54
    Aim:
        Reproduce table formatting required for src/r/deseq2_edger.R script from table produced by featureCounts tool.
    Test:
        out/awk/extract_counts_from_featureCounts/featureCounts/gtfGRCm38_texon_ggene_id_khan_2017_02_14.tsv"
        out/awk/extract_counts_from_featureCounts/subread/featureCounts_gtf-hg38-ensembl_t-exon_g-transcript_id_bam-hg38-RNA-thymus.tsv
    """
    input:
        tsv="out/{filler}.tsv"
        #tsv="out/featureCounts/gtf{gtf_descriptor}_t{feature_type}_g{attribute_type}_{bam_list}.tsv"
    output:
        tsv="out/awk/extract_counts_from_featureCounts/{filler}.tsv"
    shell:
        """
        cut -f 1,7- {input.tsv} | awk 'NR > 1' | \
            awk 'BEGIN{{FS=OFS="\t"}}; {{ for (i=2; i<=NF; i++) {{gsub(".*/","",$i)}}; print }}' | \
            awk '{{if(NR==1){{gsub(".bam","",$0); print}}else{{print}} }}' > {output.tsv} 
        """

