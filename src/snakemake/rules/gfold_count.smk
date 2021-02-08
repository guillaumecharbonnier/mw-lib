rule gfold_count_bed:
    """
    Created:
        2017-03-15 11:41:36
    Aim:
        Testing using a bed instead of a gtf as reference for counts.
    Suggested preceding rules:
        STAR, samtools_sam_to_bam
    Suggested following rules:
        gfold_diff.
    Note:
        Old export before using conda.
        #export LD_LIBRARY_PATH={WDIR}/{input.lib}
    Test:
        Testing what GFOLD does when given a bed with 'exon-like' dispersed in the file (because in this case these 'exons' are the members of a repeat family around the genome):        out/gfold/count_bed-repeat-masker/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-KO.tsv
        First test:
            out/gfold/count_bed-h4k5ac-peaks/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm10/sickle/pe_-t_sanger_-q_30/merge_lanes/run140_run141/H4K5ac-Nut-KO.tsv
    """
    input:
        bam = "out/{filler}.bam",
        bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_id])
    output:
        tsv="out/gfold/count_{bed_id}/{filler}.tsv"
    wildcard_constraints:
        bed_id="bed-[a-zA-Z0-9-]+"
    conda:
        "../envs/gfold.yaml"
    shell:
        """
        samtools view {input.bam} | \
            gfold count \
                -ann {input.bed} \
                -annf BED \
                -tag stdin \
                -o {output.tsv}
        """

rule gfold_count_gtf:
    """
    Created:
        2017-03-25 17:40:37
    Modified:
        2017-09-27 09:14:06 - input gtf function changed from input_gtf_gfold_count to input_gtf.
    Aim:
        Testing using a gtf instead of a gtf as reference for counts.
    Suggested preceding rules:
        STAR, samtools_sam_to_bam
    Suggested following rules:
        gfold_diff.
    Test:
        out/gfold/count_gtf-GRCh38-ensembl/star/pe_staridx-GRCh38-BP-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.tsv
    
    out/star/pe_staridx-GRCh38_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam· Thy3_S1·CON·PE· stranded
    out/star/pe_staridx-GRCh38_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040278.bam· Thy3_S2·CON·PE· stranded
    out/star/pe_staridx-GRCh38_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040279.bam· Thy4_S1·TRT·PE· stranded
    out/star/pe_staridx-GRCh38_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040280.bam· Thy4_S4·TRT·PE· stranded


    """
    input:
        bam = "out/{filler}.bam",
        gtf = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_id])
    output:
        tsv="out/gfold/count_{gtf_id}/{filler}.tsv"
    wildcard_constraints:
        gtf_id="gtf-[a-zA-Z0-9-]+"
    conda:
        "../envs/gfold.yaml"
    shell:
        """
        samtools view {input.bam} | \
            gfold count \
                -ann {input.gtf} \
                -annf GTF \
                -tag stdin \
                -o {output.tsv}
        """


