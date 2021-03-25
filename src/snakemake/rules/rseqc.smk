rule rseqc_geneBody_coverage:
    """
    Doc:
        http://rseqc.sourceforge.net/#genebody-coverage-py
    Test:
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
        bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id])
    output:
        pdf = "out/rseqc/geneBody_coverage_{bed_id}/{filler}.geneBodyCoverage.curves.pdf"
    params:
        outprefix = "out/rseqc/geneBody_coverage_{bed_id}/{filler}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.outprefix}"

rule rseqc_geneBody_coverage_bam_list:
    """
    Note:
        example of required bed:
        bed="out/sed/remove_chr/awk/extract_main_chr/gunzip/wget/sourceforge_rseqc/BED/Mouse_Mus_musculus/mm10.HouseKeepingGenes.bed"
    Test:
    """
    input:
        bam = lambda wildcards: eval(mwconf['ids'][wildcards.bam_list_id]),
        # Optional: Make sure bam files have their bai.
        # bai = todo_write_a_function_that_take_bam_and_add_suffix_bai,
        bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id])
    output:
        pdf="out/rseqc/geneBody_coverage_bam_list_{bed_id}/{bam_list_id}.geneBodyCoverage.curves.pdf"
    params:
        outprefix="out/rseqc/geneBody_coverage_bam_list_{bed_id}/{bam_list_id}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        BAM_LIST=`echo "{input.bam}" | tr ' ' ','`
        geneBody_coverage.py -i $BAM_LIST -r {input.bed} -o {params.outprefix}
        """

rule rseqc_infer_experiment:
    """
    Created:
        2019-02-05 23:25:41
    Doc:
        http://rseqc.sourceforge.net/#infer-experiment-py
    Test:
        out/rseqc/infer_experiment/bed-hg38-GENCODE-knownGene/ln/updir/mw/inp/bam/hg19_RNA-Seq_thymus/ISP.bam.txt
        Does not work for this bam and I do not know why:
        out/rseqc/infer_experiment/bed-hg38-GENCODE-knownGene/out/star/pe_staridx-GRCh38_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam.txt
    """
    input:
        sam_or_bam="out/{filler}.{sam_or_bam}",
        bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id])
    output:
        txt = "out/rseqc/infer_experiment/{bed_id}/{filler}.{sam_or_bam}.txt"
    log:
              "out/rseqc/infer_experiment/{bed_id}/{filler}.{sam_or_bam}.log"
    benchmark:
              "out/rseqc/infer_experiment/{bed_id}/{filler}.{sam_or_bam}.benchmark.tsv"
    wildcard_constraints:
        sam_or_bam='sam|bam',
        bed_id='bed-[a-zA-Z0-9-]+'
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.sam_or_bam} > {output.txt} 2> {log}"

