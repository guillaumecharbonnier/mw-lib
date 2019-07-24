rule snakemake_flowcharts_legacy:
    """
    Modified:
        2017-05-03 12:11:31
    Aim:
        Draw flowcharts for the given snakefile: dag and rulegraph.
    Test:

    """
    input:
        workflow="src/snakemake/workflows/{id}.snakefile"
        #workflow = "out/{id}.snakefile"
    output:
        rot_dag_dot   = "out/snakemake/rot_dag/{id}.dot",
        rot_dag_pdf   = "out/snakemake/rot_dag/{id}.pdf",
        rot_dag_png   = "out/snakemake/rot_dag/{id}.png",
        dag_dot       = "out/snakemake/dag/{id}.dot",
        dag_pdf       = "out/snakemake/dag/{id}.pdf",
        dag_png       = "out/snakemake/dag/{id}.png",
        rulegraph_dot = "out/snakemake/rulegraph/{id}.dot",
        rulegraph_pdf = "out/snakemake/rulegraph/{id}.pdf",
        rulegraph_png = "out/snakemake/rulegraph/{id}.png"
    conda:
        "../envs/snakemake.yaml"
    threads:
        1
    shell:
        """
        snakemake \
            -s {input.workflow} \
            --rulegraph > {output.rulegraph_dot}
        dot {output.rulegraph_dot} -Tpdf -o {output.rulegraph_pdf}
        dot {output.rulegraph_dot} -Tpng -o {output.rulegraph_png}

        snakemake \
            -s {input.workflow} \
            --dag > {output.dag_dot}
        dot {output.dag_dot} -Tpdf -o {output.dag_pdf}
        dot {output.dag_dot} -Tpng -o {output.dag_png}

        cat {output.dag_dot} | perl -pe 's/graph\[/graph\[rankdir='LR', /' > {output.rot_dag_dot}
        dot {output.rot_dag_dot} -Tpdf -o {output.rot_dag_pdf}
        dot {output.rot_dag_dot} -Tpng -o {output.rot_dag_png}
        """
