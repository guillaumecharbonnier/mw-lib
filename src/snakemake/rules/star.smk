rule tests_star_scipio:
    input:
        expand(
            "out/star/scipio_staridx-GRCh38-ensembl-r110_gtf-GRCh38-ensembl-r110/ln/alias/sst/all_samples/fastq/{sample}.bam",
            sample = [
                "HSC_BplusFLT3_Scipio",
                "HSC_BplusGCSF_Scipio",
                "HSC_BplusMCSF_Scipio",
                "HSC_B_Scipio",
                "MLP_BplusFLT3_Scipio",
                "MLP_B_Scipio"
            ]
        )

rule star_scipio:
    """
    Created:
        2025-02-19 16:12:59
    Aim:
        Apply STARsolo according to Scipio guidelines: Scipio-bioscience_Cytonaut-Local-UserGuide_v1.3.pdf
    Test:
        out/star/scipio_staridx-GRCh38-ensembl-r110_gtf-GRCh38-ensembl-r110/ln/alias/sst/all_samples/fastq/MLP_B_Scipio.bam
    """
    input:
        fwd = "out/{filler}_1.fastq.gz",
        rev = "out/{filler}_2.fastq.gz",
        white_list = "out/sh/generate_scipio_white_list/{filler}_1.white_list.tsv",
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        gtf   = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_id])
    output:
        bam    = "out/{tool}{extra}_{index_id}_{gtf_id}/{filler}.bam",
        log    = "out/{tool}{extra}_{index_id}_{gtf_id}/{filler}/Log.final.out"
    log:
        "out/{tool}{extra}_{index_id}_{gtf_id}/{filler}.log"
    params:
        outdir = "out/{tool}{extra}_{index_id}_{gtf_id}/{filler}",
        # Adding tmp dir here because I get FIFO error when outdir is inside the rclone mounted folder.
        tmpdir = "/tmp/{tool}{extra}_{index_id}_{gtf_id}/{filler}",
        extra = params_extra
    wildcard_constraints:
        tool="star/scipio",
        gtf_id="gtf-[a-zA-Z0-9-]+",
        star_index_id="staridx-[a-zA-Z0-9-]+"
    conda:
        "../envs/star.yaml"
    threads:
        24
    shell:
        """
        ( WDIR=`pwd`
        # Clean tmp dir and outdir if previous execution was interrupted
        rm -rf {params.tmpdir} {params.outdir}
        mkdir -p {params.outdir} `dirname {params.tmpdir}`
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{input.index} \
            --readFilesIn $WDIR/{input.rev} $WDIR/{input.fwd} \
            --sjdbGTFfile $WDIR/{input.gtf} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outFilterMultimapNmax 10 \
            --outTmpDir {params.tmpdir} \
            --twopassMode Basic --quantMode GeneCounts \
            --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --soloType CB_UMI_Simple --soloFeatures Gene --soloCBmatchWLtype Exact \
            --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
            --soloCBwhitelist $WDIR/{input.white_list} --soloCBstart 1 \
            --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 13 \
            {params.extra}

            # --soloBarcodeReadLength 31 
            # was inside the Scipio guideline but it is weird according of the other arguments.
            # Testing without it to see if it avoid critical error


            # Need to check relevance of this:
            # --outFileNamePrefix %sample_name%_ 

        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
        rm -rf {params.outdir}/_STARgenome) &> {log}
        """

rule star_pe_extra:
    """
    Modified:
        2017-03-25 14:18:31 - updated patterns for inputs
    Test:
       out/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam
       out/star/pe_fastq.gz_to_bam_--quantMode_GeneCounts_staridx-GRCm38-ensembl_gtf-GRCm38-ensembl/sickle/pe_-t_sanger_-q_30/pysradb_download_gsm_pe/GSM4475443.bam
    """
    input:
        fwd = "out/{filler}_1.{extin}",
        rev = "out/{filler}_2.{extin}",
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        gtf   = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_id])
    output:
        bam    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}.{extout}",
        log    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}/Log.final.out"
    log:
        "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}.log"
    params:
        outdir = "out/star/pe_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}",
        tmpdir = "/tmp/star/pe_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}",
        extra = params_extra,
        extin = lambda wildcards: "--readFilesCommand zcat" if wildcards.extin == 'fastq.gz' else "",
        extout = lambda wildcards: "--outSAMtype BAM SortedByCoordinate" if wildcards.extout == 'bam' else ""
    wildcard_constraints:
        tool="star/pe",
        gtf_id="gtf-[a-zA-Z0-9-]+",
        star_index_id="staridx-[a-zA-Z0-9-]+",
        extin="fastq|fastq.gz",
        extout="bam|sam"
    conda:
        "../envs/star.yaml"
    threads:
        8
    shell:
        """
        ( WDIR=`pwd`
        rm -rf {params.tmpdir} {params.outdir}
        mkdir -p {params.outdir} `dirname {params.tmpdir}`

        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{input.index} \
            --readFilesIn $WDIR/{input.fwd} $WDIR/{input.rev} \
            --sjdbGTFfile $WDIR/{input.gtf} \
            --runThreadN {threads} \
            --outTmpDir {params.tmpdir} \
            {params.extin} {params.extout} {params.extra}
        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
        rm -rf {params.outdir}/_STARgenome ) &> {log}
        """

rule star_se_extra:
    """
    Created:
        2017-02-14 11:52:08
    Aim:
        Align RNA-Seq reads with STAR.
    Doc:
        https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
    Test:
        out/star/se_fastq.gz_to_bam_standard_staridx-GRCh38_gtf-GRCh38-ensembl/sra-tools/fastq-dump_se/SRR646457.bam
        out/star/se_fastq.gz_to_bam_standard_staridx-GRCm38-ensembl_gtf-GRCm38-ensembl/sra-tools/fastq-dump_se/SRR8689837.bam
        out/star/se_fastq.gz_to_bam_standard_staridx-hg19-ensembl_gtf-hg19-ensembl/sickle/se_-t_sanger_-q_20/ln/alias/sst/RNA/run150/fastq/N2_3.bam
    """
    input:
        fwd="out/{filler}.{extin}",
        index = lambda wildcards: mwconf['ids'][wildcards.index_id],
        gtf   = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_id])
    output:
        bam    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}.{extout}",
        log    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}/Log.final.out"
    params:
        outdir = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}",
        tmpdir = "/tmp/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}",
        extra = params_extra,
        extin = lambda wildcards: "--readFilesCommand zcat" if wildcards.extin == 'fastq.gz' else "",
        extout = lambda wildcards: "--outSAMtype BAM SortedByCoordinate" if wildcards.extout == 'bam' else "",
    wildcard_constraints:
        tool="star/se",
        gtf_id="gtf-[a-zA-Z0-9-]+",
        star_index_id="staridx-[a-zA-Z0-9-]+",
        extin="fastq|fastq.gz",
        extout="bam|sam"
    conda:
        "../envs/star.yaml"
    threads:
        8
    shell:
        """
        WDIR=`pwd`
        rm -rf {params.tmpdir} {params.outdir}
        mkdir -p {params.outdir} `dirname {params.tmpdir}`
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{input.index} \
            -c --readFilesIn $WDIR/{input.fwd}\
            --sjdbGTFfile $WDIR/{input.gtf} \
            --runThreadN {threads} \
            --outTmpDir {params.tmpdir} \
            {params.extin} {params.extout} {params.extra}
        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
        rm -rf {params.outdir}/_STARgenome
        """


rule star_se_legacy:
    """
    Created:
        2017-02-14 11:52:08
    """
    input:
        fwd="out/{filler}.fastq.gz",
        index="out/star/build_index/{index}/Genome",
        gtf="out/ln/annotation_ensembl/{index}.gtf"
    output:
        bam="out/star/se_{index}/{filler}.bam",
    params:
        genomedir="out/star/build_index/{index}",
        outdir="out/star/se_{index}/{filler}"
    conda:
        "../envs/star.yaml"
    threads:
        8
    shell:
        """
        WDIR=`pwd`
        mkdir -p {params.outdir}
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{params.genomedir} \
            --readFilesCommand zcat \
            -c --readFilesIn $WDIR/{input.fwd}\
            --runThreadN {threads} \
            --sjdbGTFfile $WDIR/{input.gtf} \
            --outFilterMismatchNoverLmax 0.08 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 1 \
            --genomeLoad NoSharedMemory
        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
        rm -rf {params.outdir}/_STARgenome
        """


rule star_pe_index_legacy:
    """
    Modified:
        2017-03-25 14:18:31 - updated patterns for inputs
    Test:
       out/star/pe_GRCh38/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam
    """
    input:
        fwd="out/{filler}_1.fastq.gz",
        rev="out/{filler}_2.fastq.gz",
        index="out/star/build_index/{index}/Genome",
        gtf="out/ln/annotation_ensembl/{index}.gtf"
    output:
        bam="out/star/pe_{index}/{filler}.bam",
        log="out/star/pe_{index}/{filler}/Log.final.out"
    params:
        genomedir="out/star/build_index/{index}",
        outdir="out/star/pe_{index}/{filler}"
    conda:
        "../envs/star.yaml"
    threads:
        8
    shell:
        """
        WDIR=`pwd`
        mkdir -p {params.outdir}
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{params.genomedir} \
            --readFilesIn $WDIR/{input.fwd} $WDIR/{input.rev} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --sjdbGTFfile $WDIR/{input.gtf} \
            --outFilterMismatchNoverLmax 0.08 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 1 \
            --genomeLoad NoSharedMemory
        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
        rm -rf {params.outdir}/_STARgenome
        """

rule star_pe_index_outFilterMultimapNmax_legacy:
    """
    Modified:
        2017-03-27 18:15:38 - Forked in order to increase the parameter 'outFilterMultimapNmax' which should allow to map the ~5% of 'too many loci'. 
    Test:
       out/star/pe_GRCm38_outFilterMultimapNmax-1_outReadsUnmapped-Fastx/sickle/pe_-t_sanger_-q_20/inp/fastq/run176/RNA-P-Nut-WT-Rep1.bam.bai 
    """
    input:
        fwd="out/{filler}_1.fastq",
        rev="out/{filler}_2.fastq",
        index="out/star/build_index/{index}/Genome",
        gtf="out/ln/annotation_ensembl/{index}.gtf"
    output:
        bam="out/star/pe_{index}_outFilterMultimapNmax-{outFilterMultimapNmax}/{filler}.bam",
        unmapped_pair1="out/star/pe_{index}_outFilterMultimapNmax-{outFilterMultimapNmax}/{filler}/unmapped_pair1.fastq",
        unmapped_pair2="out/star/pe_{index}_outFilterMultimapNmax-{outFilterMultimapNmax}/{filler}/unmapped_pair2.fastq",
        log="out/star/pe_{index}_outFilterMultimapNmax-{outFilterMultimapNmax}/{filler}/Log.final.out" # Contain useful stats on uniquely mapped, too mani loci, etc...
    params:
        genomedir="out/star/build_index/{index}",
        outdir="out/star/pe_{index}_outFilterMultimapNmax-{outFilterMultimapNmax}/{filler}"
    wildcard_constraints:
        outFilterMultimapNmax="[0-9]+",
        #outReadsUnmapped="Fastx|None"
    conda:
        "../envs/star.yaml"
    threads:
        8
    shell:
        """
        WDIR=`pwd`
        mkdir -p {params.outdir}
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{params.genomedir} \
            --readFilesIn $WDIR/{input.fwd} $WDIR/{input.rev} \
            --runThreadN {threads} \
            --sjdbGTFfile $WDIR/{input.gtf} \
            --outFilterMismatchNoverLmax 0.08 \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax {wildcards.outFilterMultimapNmax} \
            --outReadsUnmapped Fastx \
            --genomeLoad NoSharedMemory

        ln -f Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
        ln -f Unmapped.out.mate1 $WDIR/{output.unmapped_pair1}
        ln -f Unmapped.out.mate2 $WDIR/{output.unmapped_pair2}
        """

rule star_build_index_tests:
    input:
        "out/star/build_index/fa-genome-GRCh38_gtf-GRCh38-merge-attr-retrieve-ensembl",
        "out/star/build_index/fa-genome-GRCh38-ensembl-r110_gtf-GRCh38-merge-attr-retrieve-ensembl-r110"
        # "out/star/build_index/fa-genome-GRCh38_gtf-GRCh38-merge-attr-retrieve-ensembl/Genome"

rule star_build_index:
    """
    Modified:
        2017-02-16 15:11:27
        2017-04-22 17:46:35 - Added tmpDir and mv Log.out.
    TODO:
        mv solution is still not perfect because log may be overwritten if multiple execution are done in parallel.
    Note:
        rm -rf tmp-dir before running because of this error:
        EXITING because of fatal ERROR: could not make temporary directory: out/star/build_index/GRCh38/tmp
        SOLUTION: (i) please check the path and writing permissions
        (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR
    """
    input:
        gtf  = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_id]),
        fa   = lambda wildcards: eval(mwconf['ids'][wildcards.fa_id])
    output:
        directory(
            "out/star/build_index/{fa_id}_{gtf_id}"
        )#,
        # "out/star/build_index/{fa_id}_{gtf_id}/Genome"
        # multiext(
        #     "out/star/build_index/{fa_id}_{gtf_id}/",
        #     "genomeParameters.txt",
        #     "Genome"
        # )
        # genomeParameters = "out/star/build_index/{fa_id}_{gtf_id}/genomeParameters.txt",
        # genome           = "out/star/build_index/{fa_id}_{gtf_id}/Genome",
    params:
        outdir           = "out/star/build_index/{fa_id}_{gtf_id}",
        tmp              = "out/star/build_index/{fa_id}_{gtf_id}_tmp"
    cache:
        "omit-software"
    conda:
        "../envs/star.yaml"
    threads:
        4
    shell:
        """
        rm -rf {params.tmp}

        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.outdir} \
            --outTmpDir {params.tmp} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf}

        # Because Star produce a log into the wdir.
		# Edit: 17/09/2020: STAR now write log in output folder. mv no more required
        #mv Log.out {params.outdir}/Log.out
        """

rule star_build_index_legacy:
    """
    Modified:
        2017-02-16 15:11:27
        2017-04-22 17:46:35 - Added tmpDir and mv Log.out.
    TODO:
        mv solution is still not perfect because log may be overwritten if multiple execution are done in parallel.
    Note:
        rm -rf tmp-dir before running because of this error:
        EXITING because of fatal ERROR: could not make temporary directory: out/star/build_index/GRCh38/tmp
        SOLUTION: (i) please check the path and writing permissions
        (ii) if you specified --outTmpDir, and this directory exists - please remove it before running STAR

    """
    input:
        gtf    = "out/ln/annotation_ensembl/{index}.gtf",
        genome = "out/cat/assembly_ensembl/{index}.fa"
    output:
        genomeParameters = "out/star/build_index/{index}/genomeParameters.txt",
        genome = "out/star/build_index/{index}/Genome",
    params:
        dir = "out/star/build_index/{index}",
        tmp = "out/star/build_index/{index}/tmp"
    conda:
        "../envs/star.yaml"
    threads:
        4
    shell:
        """
        rm -rf {params.tmp}

        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.dir} \
            --outTmpDir {params.tmp} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf}

        # Because Star produce a log into the wdir.
        #mv Log.out {params.dir}/Log.out
        """

