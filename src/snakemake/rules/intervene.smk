rule intervene_test:
    """
    Created:
        2018-02-13 13:59:48
    Aim:
        Test intervene on Thymus data.
    Test:
        out/intervene/hg38-macs2-peaks-H3K27ac-thymus-stage-exclusive/Intervene_venn.pdf
        out/intervene/hg38-macs2-peaks-H3K27ac-thymus-no-peaks-in-input/Intervene_venn.pdf
    """
    input:
        bed_list = lambda wildcards: eval(config['ids'][wildcards.bed_list_id])
    output:
        "out/intervene/{bed_list_id}/Intervene_venn.pdf",
        "out/intervene/{bed_list_id}/Intervene_upset.pdf",
    params:
        outdir="out/intervene/{bed_list_id}"
    conda:
        "../envs/intervene.yaml"
    shell:
        """
        NUMBER_OF_BED=`echo "{input.bed_list}" | wc -w`
        if [ "$NUMBER_OF_BED" -lt "6" ] ;
        then
            intervene venn -i {input.bed_list} --output {params.outdir}
        else
            echo "Too many bed files, Venn not possible"
        fi
        intervene upset -i {input.bed_list} --output {params.outdir}
        # Pairwise throws error if there are no intersection so I comment it for the moment.
        #intervene pairwise -i {input.bed_list} --output {params.outdir}
        """

