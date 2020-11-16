rule star_pe_extra:
    """
    Modified:
        2017-03-25 14:18:31 - updated patterns for inputs
    Test:
       out/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/sickle/pe_-t_sanger_-q_30/sra-tools/fastq-dump_pe/SRR2040277.bam
    """
    input:
        fwd = "out/{filler}_1.{extin}",
        rev = "out/{filler}_2.{extin}",
        index = lambda wildcards: config['ids'][wildcards.index_id],
        gtf   = lambda wildcards: eval(config['ids'][wildcards.gtf_id])
    output:
        bam    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}.{extout}",
        log    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}/Log.final.out"
    log:
                 "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}.log"
    params:
        outdir = "out/star/pe_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}",
        extra = params_extra,
        extin = lambda wildcards: "--readFilesCommand zcat" if wildcards.extin == 'fastq.gz' else "",
        extout = lambda wildcards: "--outSAMtype BAM SortedByCoordinate" if wildcards.extout == 'bam' else "",
        genomedir = lambda wildcards: os.path.dirname(config['ids'][wildcards.index_id])
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
        mkdir -p {params.outdir}
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{params.genomedir} \
            --readFilesIn $WDIR/{input.fwd} $WDIR/{input.rev} \
            --sjdbGTFfile $WDIR/{input.gtf} \
            --runThreadN {threads} \
            {params.extin} {params.extout} {params.extra}
        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam} ) &> {log}
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
        index = lambda wildcards: config['ids'][wildcards.index_id],
        gtf   = lambda wildcards: eval(config['ids'][wildcards.gtf_id])
    output:
        bam    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}.{extout}",
        log    = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}/Log.final.out"
    params:
        outdir = "out/{tool}_{extin}_to_{extout}{extra}_{index_id}_{gtf_id}/{filler}",
        extra = params_extra,
        extin = lambda wildcards: "--readFilesCommand zcat" if wildcards.extin == 'fastq.gz' else "",
        extout = lambda wildcards: "--outSAMtype BAM SortedByCoordinate" if wildcards.extout == 'bam' else "",
        genomedir = lambda wildcards: os.path.dirname(config['ids'][wildcards.index_id])
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
        mkdir -p {params.outdir}
        cd {params.outdir}
        STAR \
            --genomeDir $WDIR/{params.genomedir} \
            -c --readFilesIn $WDIR/{input.fwd}\
            --sjdbGTFfile $WDIR/{input.gtf} \
            --runThreadN {threads} \
            {params.extin} {params.extout} {params.extra}
        mv Aligned.sortedByCoord.out.bam $WDIR/{output.bam}
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
    Test:
        out/star/build_index/fa-genome-GRCh38-Blueprint_gtf-GRCh38-ensembl/Genome
    """
    input:
        gtf  = lambda wildcards: eval(config['ids'][wildcards.gtf_id]),
        fa   = lambda wildcards: eval(config['ids'][wildcards.fa_id])
    output:
        genomeParameters = "out/star/build_index/{fa_id}_{gtf_id}/genomeParameters.txt",
        genome           = "out/star/build_index/{fa_id}_{gtf_id}/Genome",
    params:
        outdir           = "out/star/build_index/{fa_id}_{gtf_id}",
        tmp              = "out/star/build_index/{fa_id}_{gtf_id}/tmp"
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

