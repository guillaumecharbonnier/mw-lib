# example of working command
# Load a recent bcl2fastq env before running
# conda activate bcl2fastq
#/gpfs/tgml/apps/cellranger-3.1.0/cellranger mkfastq --run=/gpfs/tgml/reads/bcl/Run_243_190213_NS500637_0150_AHNJM5BGX7 --id=Run_243 --csv=/gpfs/tgml/reads/bcl/Run_243_190213_NS500637_0150_AHNJM5BGX7/SampleSheet.csv --output-dir=/gpfs/tgml/nin/test_cellranger
###
# I encountered some trouble with the param extra. So for now, the localcores and localmem will be set manually in the rule. 
###
#
# 2021-03-06 Update cellranger to v 6.0.0


rule cellranger_mkfastq:
    """
    Doc:
        https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq 
    Test:
        out/cellranger/mkfastq/gpfs/tgml/reads/bcl/Run_309_200224_NS500637_0216_AH2J22BGXF/SampleSheet.csv
    """
    input:
        xml="/{filler}/RunInfo.xml",
        csv="out/{tool}{extra}/{filler}/SampleSheet.csv"
    output:
        xml="out/{tool}{extra}/{filler}/RunInfo.xml",
        html="out/{tool}{extra}/{filler}/Reports/html/tree.html",
    conda:
        "../envs/bcl2fastq.yml"
    params:
        extra = params_extra
    threads:
        8
    wildcard_constraints:
        tool="cellranger/mkfastq"
    log:
        mkfastq="out/{tool}{extra}/{filler}/mkfastq_log"
    shell:'''
        (
        # 2021-05-27: Useless condition, snakemake don't check this before running the rule
        #Condition to not copy at each execution
        #if [[ -f {output.xml} ]]; then
        #    echo "RunInfo.xml exists."
        #    break
        #else
        #    cp {input.xml} {output.xml}      
        #fi
      
        INDIR=`dirname {input.xml}`
        OUTDIR=`dirname {input.csv}`
        EXP=`grep Experiment {input.csv} | cut -d "," -f 2`
        RUN=`echo {input.csv} | awk -F"/" '{{print $8}}' | cut -d "_" -f 1,2`
        # Trick to get the relative path to cellranger from output-dir
        CELLRANGER_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('../apps/cellranger-6.0.0', '${{OUTDIR}}'))"`
        export PATH=$CELLRANGER_RELATIVE_PATH_TO_OUTPUT:$PATH
        # Same trick as before to get relative path to fastq file from output-dir
        INDIR_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('${{INDIR}}', '${{OUTDIR}}'))"`
        cd ${{OUTDIR}}
        #alias bcl2fastq="/gpfs/tgml/mw-sst/.snakemake/conda/3b33e592/bin/bcl2fastq"
        cellranger mkfastq {params.extra} --run=${{INDIR_RELATIVE_PATH_TO_OUTPUT}} --id=${{RUN}} --csv=./SampleSheet.csv
        cp -r ${{RUN}}/outs/fastq_path/Reports .
        ) &> {log}
        # go back to mw-sst folder
        #cd /gpfs/tgml/mw-sst
        '''


# Created on 02/09/2020
# Generate count for fastq files generated by the cellranger_mkfastq rule above.
rule cellranger_count:
    """
    Doc:
        https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
    Test:
        out/tar/xvzf/wget/https/cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A/done
    Output:
        out/cellranger/count_mm10-2020-A/cellranger/mkfastq/gpfs/tgml/reads/bcl/Run_243_190213_NS500637_0150_AHNJM5BGX7/web_summary.html
        out/cellranger/count_GRCh38-2020-A/cellranger/mkfastq/gpfs/tgml/reads/bcl/Run_309_200224_NS500637_0216_AH2J22BGXF/web_summary.html
    """
    input:
        ref="out/tar/xvzf_genome_cellranger/wget/https/cf.10xgenomics.com/supp/cell-exp/refdata-gex-{assembly}/done",
        xml="out/{filler}/RunInfo.xml"
    output:
        html="out/{tool}{extra}_{assembly}/{filler}/process_done"
    params:
        extra = params_extra
    threads:
        16
    wildcard_constraints:
        tool="cellranger/count",
        assembly="[A-Za-z0-9-]+"
    log:
        "out/{tool}{extra}_{assembly}/{filler}/log"
    shell:'''
        (
        INDIR=`dirname {input.xml}`
        OUTDIR=`dirname {output.html}`
        REF=`dirname {input.ref}`
        SAMPLESHEET=${{INDIR}}"/SampleSheet.csv"
        EXP=`grep Experiment ${{SAMPLESHEET}} | cut -d , -f 2`
        # Extract sample names from samplesheet
        SAMPLE=`grep S0 ${{SAMPLESHEET}} | grep -v -E "HTO|ADT" | cut -d , -f 2`
        RUN=`grep Experiment ${{SAMPLESHEET}} | cut -d , -f 2 | cut -d _ -f 1,2`
        # Trick to get the relative path to cellranger from output-dir
        CELLRANGER_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('../apps/cellranger-6.0.0', '${{OUTDIR}}'))"`
        export PATH=$CELLRANGER_RELATIVE_PATH_TO_OUTPUT:$PATH
        # Same trick as before to get relative path to fastq file from output-dir
        INDIR_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('${{INDIR}}', '${{OUTDIR}}'))"`
        # And also to have relative path to ref from output-dir
        REF_RELATIVE_PATH_TO_OUTPUT=`python -c "import os.path; print(os.path.relpath('${{REF}}', '${{OUTDIR}}'))"`
        # Move into the outdir to have cellranger count output in the correct folder 
        cd ${{OUTDIR}}
        #cellranger count {params.extra} --id=${{RUN}} --fastqs=${{INDIR_RELATIVE_PATH_TO_OUTPUT}} --transcriptome=${{REF_RELATIVE_PATH_TO_OUTPUT}} --project=${{EXP}} --sample=${{SAMPLE}}
        #cellranger count {params.extra} --id=${{RUN}} --fastqs=${{INDIR_RELATIVE_PATH_TO_OUTPUT}}/${{RUN}}/outs/fastq_path --transcriptome=${{REF_RELATIVE_PATH_TO_OUTPUT}} --sample=${{SAMPLE}}
        cellranger count {params.extra} --id=${{RUN}} --fastqs=${{INDIR_RELATIVE_PATH_TO_OUTPUT}}/${{RUN}}/outs/fastq_path --transcriptome=${{REF_RELATIVE_PATH_TO_OUTPUT}} 
        touch process_done 
        )&> {log}
        # go back to mw-sst folder
        cd /gpfs/tgml/mw-sst
    '''
