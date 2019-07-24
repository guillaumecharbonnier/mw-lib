rule r_featureCounts_to_bed_of_transcript_classes:
    """
    Created:
        2017-09-20 15:45:54
    Test:
        out/r/featureCounts_to_bed_of_transcript_classes/featureCounts/gtf-longs-protein-coding-merge-attr-GRCm38_t-exon_g-merge_gene_id_name_bam-h2al2-rnaseq/done
        out/r/featureCounts_to_bed_of_transcript_classes_pattern-RNA-C-H2AL2-WT/featureCounts/gtf-longs-protein-coding-merge-attr-GRCm38_t-exon_g-merge_gene_id_name_bam-h2al2-rnaseq/done
    """
    input:
        rscript="opt/miniconda/envs/deseq2_edger/bin/Rscript", # I just need getopt for this script
        script="src/r/draft/featureCounts_to_bed_of_transcript_classes.R",
        tsv="out/{filler}.tsv"
    output:
        done="out/r/featureCounts_to_bed_of_transcript_classes_pattern-{pattern}/{filler}/done",
        bed=expand(
            "out/r/featureCounts_to_bed_of_transcript_classes_pattern-{{pattern}}/{{filler}}/{classes}.bed",
            classes=[
                "RPKM_mean_gtq0.9",
                "RPKM_mean_gtq0.8_ltq0.9",
                "RPKM_mean_gtq0.7_ltq0.8",
                "RPKM_mean_gtq0.6_ltq0.7",
                "RPKM_mean_gtq0.5_ltq0.6",
                "RPKM_mean_gtq0.4_ltq0.5",
                "RPKM_mean_gtq0.3_ltq0.4",
                "RPKM_mean_gtq0.2_ltq0.3",
                "RPKM_mean_gtq0.1_ltq0.2",
                "RPKM_mean_ltq0.1",
                "RPKM_mean_gtQ3",
                "RPKM_mean_gtQ2_ltQ3",
                "RPKM_mean_gtQ1_ltQ2",
                "RPKM_mean_ltQ1",
                "RPKM_mean_0"])
    params:
        outdir="out/r/featureCounts_to_bed_of_transcript_classes_pattern-{pattern}/{filler}"
    wildcard_constraints:
        pattern="[a-zA-Z0-9-_]+"
    shell:
        """
        {input.rscript} {input.script} \
            --input {input.tsv} \
            --outdir {params.outdir} \
            --pattern {wildcards.pattern}


        touch {output.done}
        """

