# example of working command
# Load a recent bcl2fastq env before runnin
# conda activate bcl2fastq
#/gpfs/tgml/apps/cellranger-3.1.0/cellranger mkfastq --run=/gpfs/tgml/reads/bcl/Run_243_190213_NS500637_0150_AHNJM5BGX7 --id=Run_243 --csv=/gpfs/tgml/reads/bcl/Run_243_190213_NS500637_0150_AHNJM5BGX7/SampleSheet.csv --output-dir=/gpfs/tgml/nin/test_cellranger


rule cellranger_mkfastq:
	input:
		xml="/{filler}/RunInfo.xml",
		csv="/{filler}/SampleSheet.csv"
	output:
		xml="out/{tool}{extra}/{filler}/RunInfo.xml",
		csv="out/{tool}{extra}/{filler}/SampleSheet.csv",
		html="out/{tool}{extra}/{filler}/Reports/html/tree.html"
	conda:
		"../envs/bcl2fastq.yml"
	params:
		extra = params_extra
	threads:
		8
	wildcard_constraints:
		tool="cellranger/"
	log:
		"out/{tool}{extra}/{filler}/log"
	shell:'''
		cp {input.csv} {output.csv}
		cp {input.xml} {output.xml}
		INDIR=`dirname {input.csv}`
		OUTDIR=`dirname {output.csv}`
		RUN=`echo {input.csv} | awk -F"/" '{{print $6}}' | cut -d "_" -f 1,2`
		export PATH=../apps/cellranger-3.1.0:$PATH
		cellranger {params.extra} --run=${{INDIR}} --id=${{RUN}} --csv={input.csv} --output-dir=${{OUTDIR}} &> {log}
		mv ${{RUN}} ${{OUTDIR}}/
		mv "__"${{RUN}}.mro ${{OUTDIR}}/
		'''
