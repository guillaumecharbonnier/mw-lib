rule r_greatr_yaml_new:
    """
    Created:
        2019-01-08 11:42:16
    Aim:
        Testing greatr with only a yaml containing the whole conf
    Note:
        More than one yaml exists for this bedlist:
        ['../mw-thymus/src/greatr/hg19-thymus-mca-cor-dims-1-2-mixed-thresholds.yaml', '../mw-thymus/src/greatr/hg19-thymus-mca-cor-dims-1-2.yaml']

    Test:
        out/r/greatr_bed-hg19-cpg-meth-call-clusters/ln/updir/mw-thymus/src/greatr/hg19-cpg-meth-call-clusters/done
        out/r/new_greatr_bed-hg19-thymus-mca-cor-dims-1-2-mixed-thresholds/ln/updir/mw-thymus/src/greatr/hg19-thymus-mca-cor-dims-1-2-mixed-thresholds/done
    """
    input:
        yaml = "out/{filler}.yaml",
        bed_list = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id])
    output:
        done=touch("out/r/new_greatr_{bed_list_id}/{filler}/done")
    conda:
        "../envs/r_greatr.yaml"
    shell:
        "greatr -y {input.yaml}"

rule r_greatr_yaml_old:
    """
    Created:
        2019-01-08 11:42:16
    Aim:
        Testing greatr with only a yaml containing the whole conf
    Note:
        More than one yaml exists for this bedlist:
        ['../mw-thymus/src/greatr/hg19-thymus-mca-cor-dims-1-2-mixed-thresholds.yaml', '../mw-thymus/src/greatr/hg19-thymus-mca-cor-dims-1-2.yaml']

    Test:
        out/r/greatr_yaml-hg19-active-enhancers-thymus/done
        out/r/greatr_yaml-hg19-thymus-mca-cor-dims-1-2/done out/r/greatr_yaml-hg19-thymus-mca-cor-dims-1-2-mixed-thresholds/done
    """
    input:
        #yaml='src/greatr/{bed_list_id}.yaml',
        yaml=input_yaml_greatr,
        bed_list = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id])
    output:
        done=touch("out/r/greatr_yaml-{bed_list_id}/done")
    conda:
        "../envs/r_greatr.yaml"
    shell:
        "greatr -y {input.yaml}"


#GREATR_OUTPUTS=[x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_greatr.txt","r")]
GREATR_OUTPUTS=[x.strip() for x in open("../mw-lib/src/snakemake/lists/outputs_greatr_2019-12-30_pdf_subset.txt","r")]

rule r_greatr_extra:
    """
    Created:
        2018-11-08 01:42:52
    Test:
        Obsolete (no -g arg anymore):
            out/r/great_heatmap_bed-hg19-cpg-meth-call-clusters/enrichment_tables.Rdata
            out/r/great_heatmap_g-GOBP_r-NULL_bed-hg19-h3k27ac-on-low-cpg-meth-call-clusters/enrichment_tables.Rdata
            out/r/great_heatmap__g-GOBP_m-bed-hg19-alex-test/enrichment_tables.done
            out/r/great_heatmap__g-GOBP__m-Hyper_Fold_Enrichment,Binom_Bonf_PValue,Hyper_Bonf_PValue,Post_Filter_Binom_Rank__s-greater,lower,lower,lower__t-2,0.05,0.05,5_bed-hg19-alex-test/enrichment_tables.done
            out/r/greatr_-g_MSDBP_bed-hg19-cpg-meth-call-clusters/enrichment_tables.done
            out/r/greatr_-g_GOBP_bed-hg19-cpg-meth-call-clusters/enrichment_tables.done
            out/r/greatr_-g_GOBP_-a_mm10_bed-mm10-ChromHMM-test19/enrichment_tables.done
            out/r/greatr_-g_GOBP_-a_hg19_bed-hg19-encode-broad-h3k4me3/enrichment_tables.done
        Current:
            out/r/greatr_bed-hg19-h3k27ac-on-all-cpg-hypometh-call-clusters/done
            out/r/greatr_bed-hg19-test-coord-cor-with-mca-dims/done
            out/r/greatr_bed-hg19-atac-thymus/done
        Note: I moved the "bed-" prefix from output to wildcard bed_list_id so examples above may only work in bed_list_id is modified in accordance.
            out/r/greatr_bed-hg19-active-enhancers-thymopoiesis-tall-samples/done
            out/r/greatr_bed-hg19-active-enhancers-interesting-classes/done
    """
    input:
        bed_list = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id])
    output:
        #tables="out/r/great_heatmap{extra}_bed-{bed_list_id}/enrichment_tables.Rdata"
        done=touch("out/{tool}{extra}_{bed_list_id}/done"),
        pdf=touch(expand("out/{{tool}}{{extra}}_{{bed_list_id}}/{outfile}", outfile=GREATR_OUTPUTS))
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "r/greatr",
        #extra = "[^yaml].*"
    conda:
        "../envs/r_greatr.yaml"
    shell:
        """
        FILES=`echo {input.bed_list} | sed 's/ /,/g'`
        OUTDIR=`dirname {output.done}`
        echo "OUTDIR: $OUTDIR"
        greatr -f $FILES -o $OUTDIR {params.extra}
        """
