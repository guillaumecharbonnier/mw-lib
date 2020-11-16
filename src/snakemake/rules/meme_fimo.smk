rule meme_fimo:
    """
    Created:
        2020-08-14 18:54:31
    Doc:
        http://meme-suite.org/doc/fimo.html
    Test:
    meme-HOCOMOCOv11-full-HUMAN out/tar_meme/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme

        out/meme/fimo_meme-HOCOMOCOv11-full-HUMAN/python/vcf_to_flanking_sequence_fa-genome-GRCh38-ensembl-r100/grep/extract-indel-from-vcf/bcftools/mpileup_fa-genome-GRCh38-ensembl-r100/samtools/index/samtools/sort/samtools/view_sam_to_bam/bowtie2/se_fasta_--rfg_1,1_-k_1_-q_-f_GRCh38-ensembl-r100/spades/se/ln/alias/sst/all_samples/fastq/802_BT11_H3K27ac/contigs/fimo.tsv
    """
    input:
        fasta="out/{filler}.fasta",
        meme = lambda wildcards: config['ids'][wildcards.meme_id]
    output:
        expand("out/meme/fimo_{{meme_id}}/{{filler}}/{file}", file=['fimo.html', 'fimo.tsv'])
    log:
               "out/meme/fimo_{meme_id}/{filler}/log"
    benchmark:
               "out/meme/fimo_{meme_id}/{filler}/benchmark.tsv"
    params:
        outdir="out/meme/fimo_{meme_id}/{filler}"
    wildcard_constraints:
        meme_id = "[a-zA-Z0-9-]+"
    conda:
        "../envs/meme.yaml"
    shell:
        """
        fimo --oc {params.outdir} {input.meme} {input.fasta} &> {log}
        """
