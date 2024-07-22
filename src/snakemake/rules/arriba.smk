rule arriba_extra:
    """
    Created:
        2024-03-29 20:01:56
    Note:
        The run_arriba.sh script is available here:
        https://github.com/suhrig/arriba/blob/master/run_arriba.sh
        and I could improve my snakemake workflow splitting it into the STAR rule, the arriba call and the samtools index.
    Test:
        out/arriba/pe_fq_GRCh38viral_ENSEMBL104/agent/trim_-v2/ln/alias/sst/all_samples/fastq/113281_PICCL/fusions.tsv
    """
    input:
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz",
        gtf = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_gtf_id}.gtf",
        fa = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_fa_id}.fa"
    output:
        tsv = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/fusions.tsv",
        bam = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/Aligned.sortedByCoord.out.bam",
        bai = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/Aligned.sortedByCoord.out.bam.bai"
    params:
        outdir = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}"
    wildcard_constraints:
        tool="arriba/pe_fq",
        arriba_fa_id = "GRCh38viral",
        arriba_gtf_id = "ENSEMBL104"
    conda:
        "../envs/arriba.yaml"
    threads:
        8
    shell:
        """
        ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

        # Using arriba_fa_id we can decide which suffix to use for blacklist, known_fusions and protein_domains:

        if [[ {wildcards.arriba_fa_id} == *"GRCh38"* ]] || [[ {wildcards.arriba_fa_id} == *"hg38"* ]]; then
            STEM="_hg38_GRCh38_v2.4.0"
        elif [[ {wildcards.arriba_fa_id} == *"GRCh37"* ]] || [[ {wildcards.arriba_fa_id} == *"hg19"* ]]; then
            STEM="_hg19_hs37d5_GRCh37_v2.4.0"
        else
            echo "Invalid arriba_fa_id"
            exit 1
        fi

        cd {params.outdir}

        run_arriba.sh \
            {WDIR}/out/arriba/download_references_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id}/STAR_index_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id} \
            {WDIR}/{input.gtf} \
            {WDIR}/{input.fa} \
            $ARRIBA_FILES/blacklist$STEM.tsv.gz \
            $ARRIBA_FILES/known_fusions$STEM.tsv.gz \
            $ARRIBA_FILES/protein_domains$STEM.gff3 \
            {threads} \
            {WDIR}/{input.fq1} \
            {WDIR}/{input.fq2}
        
        """

rule arriba_draw_fusions:
    """
    Created:
        2024-07-22 06:28:45
    Doc:
        https://arriba.readthedocs.io/en/latest/visualization/
    Note:
        Currently only developped to work with conda:
        The database files (protein domain track, cytobands) are located in $CONDA_PREFIX/var/lib/arriba.
    Test:
        out/arriba/pe_fq_GRCh38viral_ENSEMBL104/agent/trim_-v2/ln/alias/sst/all_samples/fastq/1396_RNA_XT_HS2/fusions.tsv
    """
    input:
        tsv = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/fusions.tsv",
        bam = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/Aligned.sortedByCoord.out.bam",
        gtf = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_gtf_id}.gtf"
    output:
        pdf = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/fusions.pdf"
    params:
        outdir = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}"
    wildcard_constraints:
        tool="arriba/pe_fq",
        arriba_fa_id = "GRCh38viral",
        arriba_gtf_id = "ENSEMBL104"
    conda:
        "../envs/arriba.yaml"
    shell:
        """
        ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

        # Using arriba_fa_id we can decide which suffix to use for cytobands and protein_domains:

        if [[ {wildcards.arriba_fa_id} == *"GRCh38"* ]] || [[ {wildcards.arriba_fa_id} == *"hg38"* ]]; then
            STEM="_hg38_GRCh38_v2.4.0"
        elif [[ {wildcards.arriba_fa_id} == *"GRCh37"* ]] || [[ {wildcards.arriba_fa_id} == *"hg19"* ]]; then
            STEM="_hg19_hs37d5_GRCh37_v2.4.0"
        else
            echo "Invalid arriba_fa_id"
            exit 1
        fi

        cd {params.outdir}

        draw_fusions.R \
            --fusions={WDIR}/{input.tsv} \
            --alignments={WDIR}/{input.bam} \
            --output={WDIR}/{output.pdf} \
            --annotation={WDIR}/{input.gtf} \
            --cytobands=$ARRIBA_FILES/cytobands$STEM.tsv \
            --proteinDomains=$ARRIBA_FILES/protein_domains$STEM.gff3
        """



# arriba \  
# -x {input.bam \  
# -g {input.gtf} \  
# -a {input.fa} \  
# -o {output.tsv} \  
# -O {output.discarded} \  
# -f blacklist

rule arriba_download_references:
    output:
        gtf = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_gtf_id}.gtf",
        fa = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_fa_id}.fa"
    params:
        outdir = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}"
    conda:
        "../envs/arriba.yaml"
    shell:
        """
        cd {params.outdir}
        $CONDA_PREFIX/var/lib/arriba/download_references.sh {wildcards.arriba_fa_id}+{wildcards.arriba_gtf_id}
        """
