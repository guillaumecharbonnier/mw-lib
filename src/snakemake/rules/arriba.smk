rule arriba_extra:
    """
    Created:
        2024-03-29 20:01:56
    Modified:
        2025-10-28 - Replaced run_arriba.sh with direct STAR/arriba calls to support --outTmpDir
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
        outdir = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}",
        tmpdir = "/tmp/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}"
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
        WDIR=`pwd`
        ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

        # Clean tmp dir if previous execution was interrupted
        rm -rf {params.tmpdir}
        mkdir -p `dirname {params.tmpdir}`

        # Using arriba_fa_id we can decide which suffix to use for blacklist, known_fusions and protein_domains:
        if [[ {wildcards.arriba_fa_id} == *"GRCh38"* ]] || [[ {wildcards.arriba_fa_id} == *"hg38"* ]]; then
            BASE="_hg38_GRCh38"
        elif [[ {wildcards.arriba_fa_id} == *"GRCh37"* ]] || [[ {wildcards.arriba_fa_id} == *"hg19"* ]]; then
            BASE="_hg19_hs37d5_GRCh37"
        else
            echo "Invalid arriba_fa_id"
            exit 1
        fi
        
        # Auto-detect the version from available files in ARRIBA_FILES
        BLACKLIST_FILE=$(ls $ARRIBA_FILES/blacklist${{BASE}}_v*.tsv.gz 2>/dev/null | head -n1)
        if [[ -z "$BLACKLIST_FILE" ]]; then
            echo "Error: No blacklist file found matching pattern: $ARRIBA_FILES/blacklist${{BASE}}_v*.tsv.gz"
            exit 1
        fi
        
        # Extract the full stem (including version) from the detected file
        STEM=$(basename "$BLACKLIST_FILE" .tsv.gz | sed "s/blacklist//")

        cd {params.outdir}

        # Run STAR with arriba settings and output to tmpdir for FIFO operations
        STAR \
            --runThreadN {threads} \
            --genomeDir $WDIR/out/arriba/download_references_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id}/STAR_index_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id} \
            --genomeLoad NoSharedMemory \
            --readFilesIn $WDIR/{input.fq1} $WDIR/{input.fq2} \
            --readFilesCommand zcat \
            --outStd BAM_Unsorted \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outBAMcompression 0 \
            --outTmpDir {params.tmpdir} \
            --outFilterMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 | \
        tee Aligned.out.bam | \
        arriba \
            -x /dev/stdin \
            -o fusions.tsv \
            -O fusions.discarded.tsv \
            -a $WDIR/{input.fa} \
            -g $WDIR/{input.gtf} \
            -b $ARRIBA_FILES/blacklist$STEM.tsv.gz \
            -k $ARRIBA_FILES/known_fusions$STEM.tsv.gz \
            -t $ARRIBA_FILES/known_fusions$STEM.tsv.gz \
            -p $ARRIBA_FILES/protein_domains$STEM.gff3

        # Sort and index BAM
        samtools sort -@ {threads} -m $((40000/{threads}))M -T tmp -O bam Aligned.out.bam > Aligned.sortedByCoord.out.bam
        rm -f Aligned.out.bam
        samtools index Aligned.sortedByCoord.out.bam
        """

# Commenting original rule for reference.
# I had to replicate run_arriba.sh content here to be able to use --outTmpDir in STAR call to avoid FIFO issue
# rule arriba_extra:
#     """
#     Created:
#         2024-03-29 20:01:56
#     Note:
#         The run_arriba.sh script is available here:
#         https://github.com/suhrig/arriba/blob/master/run_arriba.sh
#         and I could improve my snakemake workflow splitting it into the STAR rule, the arriba call and the samtools index.
#     Test:
#         out/arriba/pe_fq_GRCh38viral_ENSEMBL104/agent/trim_-v2/ln/alias/sst/all_samples/fastq/113281_PICCL/fusions.tsv
#     """
#     input:
#         fq1 = "out/{filler}_1.fastq.gz",
#         fq2 = "out/{filler}_2.fastq.gz",
#         gtf = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_gtf_id}.gtf",
#         fa = "out/arriba/download_references_{arriba_fa_id}_{arriba_gtf_id}/{arriba_fa_id}.fa"
#     output:
#         tsv = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/fusions.tsv",
#         bam = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/Aligned.sortedByCoord.out.bam",
#         bai = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}/Aligned.sortedByCoord.out.bam.bai"
#     params:
#         outdir = "out/{tool}_{arriba_fa_id}_{arriba_gtf_id}/{filler}"
#     wildcard_constraints:
#         tool="arriba/pe_fq",
#         arriba_fa_id = "GRCh38viral",
#         arriba_gtf_id = "ENSEMBL104"
#     conda:
#         "../envs/arriba.yaml"
#     threads:
#         8
#     shell:
#         """
#         ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba

#         # Using arriba_fa_id we can decide which suffix to use for blacklist, known_fusions and protein_domains:

#         if [[ {wildcards.arriba_fa_id} == *"GRCh38"* ]] || [[ {wildcards.arriba_fa_id} == *"hg38"* ]]; then
#             STEM="_hg38_GRCh38_v2.5.0"
#         elif [[ {wildcards.arriba_fa_id} == *"GRCh37"* ]] || [[ {wildcards.arriba_fa_id} == *"hg19"* ]]; then
#             STEM="_hg19_hs37d5_GRCh37_v2.5.0"
#         else
#             echo "Invalid arriba_fa_id"
#             exit 1
#         fi

#         cd {params.outdir}

#         run_arriba.sh \
#             {WDIR}/out/arriba/download_references_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id}/STAR_index_{wildcards.arriba_fa_id}_{wildcards.arriba_gtf_id} \
#             {WDIR}/{input.gtf} \
#             {WDIR}/{input.fa} \
#             $ARRIBA_FILES/blacklist$STEM.tsv.gz \
#             $ARRIBA_FILES/known_fusions$STEM.tsv.gz \
#             $ARRIBA_FILES/protein_domains$STEM.gff3 \
#             {threads} \
#             {WDIR}/{input.fq1} \
#             {WDIR}/{input.fq2}
        
#         """

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
            STEM="_hg38_GRCh38_v2.5.0"
        elif [[ {wildcards.arriba_fa_id} == *"GRCh37"* ]] || [[ {wildcards.arriba_fa_id} == *"hg19"* ]]; then
            STEM="_hg19_hs37d5_GRCh37_v2.5.0"
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
