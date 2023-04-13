# Created on 23/02/2023
# Generate fastq files using cellranger_atac_mkfastq.
rule cellranger_atac_mkfastq:
    """
    Doc:
        https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq 
    Test:
        out/cellranger-atac/mkfastq/gpfs/tgml/reads/bcl/Run_309_200224_NS500637_0216_AH2J22BGXF/SampleSheet.csv
    """
    input:
        xml="/{filler}/RunInfo.xml",
        csv="out/{tool}{extra}/{filler}/SampleSheet.csv"
    output:
        xml="out/{tool}{extra}/{filler}/RunInfo.xml",
        html="out/{tool}{extra}/{filler}/Reports/html/tree.html"
    conda:
        "../envs/bcl2fastq.yml"
    params:
        extra = params_extra
    threads:
        8
    wildcard_constraints:
        tool="cellranger-atac/mkfastq"
    log:
        mkfastq="out/{tool}{extra}/{filler}/mkfastq_log"
    shell:'''
        (
        cp {input.xml} {output.xml}
        INDIR=`dirname {input.xml}`
        OUTDIR=`dirname {input.csv}`
        EXP=`grep 'Experiment Name' {input.csv} | cut -d "," -f 2`
        RUN=`grep 'Experiment Name' {input.csv} | cut -d "," -f 2 | cut -d "_" -f 1,2`
        # Trick to get the relative path to cellranger from output-dir
        CELLRANGER_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('../apps/cellranger-atac-2.1.0', '${{OUTDIR}}'))"`
        export PATH=$CELLRANGER_RELATIVE_PATH_TO_OUTPUT:$PATH
        # Same trick as before to get relative path to fastq file from output-dir
        INDIR_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('${{INDIR}}', '${{OUTDIR}}'))"`
        cd ${{OUTDIR}}
        cellranger-atac mkfastq {params.extra} --run=${{INDIR_RELATIVE_PATH_TO_OUTPUT}} --id=${{RUN}} --csv=./SampleSheet.csv
        cp -r ${{RUN}}/outs/fastq_path/Reports .
        ) &> {log}
        '''


# Created on 23/02/2023
# Generate count for fastq files generated by the cellranger_atac_mkfastq rule above.
rule cellranger_atac_count:
    """
    Doc:
        https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
    Test:
        out/tar/xvzf/wget/https/cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A/done
    Output:
        out/cellranger_atac/count_mm10-2020-A/cellranger/mkfastq/gpfs/tgml/reads/bcl/Run_243_190213_NS500637_0150_AHNJM5BGX7/web_summary.html
        out/cellranger_atac/count_GRCh38-2020-A/cellranger/mkfastq/gpfs/tgml/reads/bcl/Run_309_200224_NS500637_0216_AH2J22BGXF/web_summary.html
    """
    input:
        ref="out/tar/xvzf_genome_cellranger/wget/https/cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-{assembly}/done",
        xml="out/{filler}/RunInfo.xml",
        csv="out/{filler}/SampleSheet.csv"
    output:
        done="out/{tool}{extra}_{assembly}/{filler}/process_done"
    params:
        extra = params_extra
    threads:
        8
    wildcard_constraints:
        tool="cellranger-atac/count",
        assembly="[A-Za-z0-9-]+"
    log:
        "out/{tool}{extra}_{assembly}/{filler}/log"
    shell:
        """
        (
        INDIR=`dirname {input.xml}`
        OUTDIR=`dirname {output.done}`
        REF=`dirname {input.ref}`
        EXP=`grep 'Experiment Name' {input.csv} | cut -d "," -f 2`
        RUN=`grep 'Experiment Name' {input.csv} | cut -d "," -f 2 | cut -d "_" -f 1,2`
        # Extract sample names from samplesheet
        # Trick to get the relative path to cellranger from output-dir
        CELLRANGER_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('../apps/cellranger-atac-2.1.0', '${{OUTDIR}}'))"`
        export PATH=$CELLRANGER_RELATIVE_PATH_TO_OUTPUT:$PATH
        # Same trick as before to get relative path to fastq file from output-dir
        INDIR_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('${{INDIR}}', '${{OUTDIR}}'))"`
        # And also to have relative path to ref from output-dir
        REF_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('${{REF}}', '${{OUTDIR}}'))"`
        # Move into the outdir to have cellranger count output in the correct folder 
        cd ${{OUTDIR}}
        cellranger-atac count {params.extra} --id=${{RUN}} --fastqs=${{INDIR_RELATIVE_PATH_TO_OUTPUT}}/${{RUN}}/outs/fastq_path --transcriptome=${{REF_RELATIVE_PATH_TO_OUTPUT}} 
        touch process_done 
        )&> {log}
        """