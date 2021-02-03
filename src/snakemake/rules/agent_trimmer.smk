rule agent_trim_pe:
    """
    Doc:
        https://www.agilent.com/cs/library/software/public/AGeNT%20ReadMe.pdf
        https://www.agilent.com/cs/library/software/public/AGeNTBestPractices.pdf
        
    Notes:
        1.  -out_loc `dirname {output.fq1}`
            does not work properly as output files are put inside input directory instead...
        2.  The "input directory" follow symlinks. As such, we should copy (and not symlink) input files into output directory before running trimmer.
        3.  The trimmer takes all available threads by default, although the RAM requirement is limited.
    Tests:
        out/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60_2.fastq.gz
        
        _L001_R1_001.fastq.gz
        MOLT4_S60_L001_R2_001.fastq.gz
    """
    input:
        fq1l1 = "out/{filler}_L001_R1_001.fastq.gz",
        fq1l2 = "out/{filler}_L002_R1_001.fastq.gz",
        fq2l1 = "out/{filler}_L001_R2_001.fastq.gz",
        fq2l2 = "out/{filler}_L002_R2_001.fastq.gz",
        agent = "out/agent/agent/agent.sh"
    output:
        fq1 = "out/{tool}{extra}/{filler}_1.fastq.gz",
        fq2 = "out/{tool}{extra}/{filler}_2.fastq.gz",
        mbc = "out/{tool}{extra}/{filler}.txt.gz",
        properties = "out/{tool}{extra}/{filler}.properties",
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "agent/trim"
    shell:
        """
        OUTDIR=`dirname {output.fq1}`
        FQ1L1=$OUTDIR/`basename {input.fq1l1}`
        FQ1L2=$OUTDIR/`basename {input.fq1l2}`
        FQ2L1=$OUTDIR/`basename {input.fq2l1}`
        FQ2L2=$OUTDIR/`basename {input.fq2l2}`

        # Input fastq symlink are dereferenced in the output directory
        # else agent will put output files into symlink targets.
        cp -rL {input.fq1l1} $FQ1L1
        cp -rL {input.fq1l2} $FQ1L2
        cp -rL {input.fq2l1} $FQ2L1
        cp -rL {input.fq2l2} $FQ2L2

        OUTPREFIX=out/{wildcards.tool}{wildcards.extra}/{wildcards.filler}

        # Output files from previous run should be removed
        # as the regex `mv` below should match only one file each.
        rm -f \
            ${{OUTPREFIX}}_L001_R1_001.fastq*_Cut_0.fastq.gz \
            ${{OUTPREFIX}}_L001_R2_001.fastq*_Cut_0.fastq.gz \
            ${{OUTPREFIX}}_L001_RN_001.fastq*_MBC_0.txt.gz \
            ${{OUTPREFIX}}_L001_RN_001.fastq*_STATS_0.properties

        {input.agent} trim -fq1 $FQ1L1,$FQ1L2 -fq2 $FQ2L1,$FQ2L2 {params.extra}

        rm -f $FQ1L1 $FQ1L2 $FQ2L1 $FQ2L2

        # Timestamped output files from agent are simplified
        mv ${{OUTPREFIX}}_L001_R2_001.fastq*_Cut_0.fastq.gz {output.fq2}
        mv ${{OUTPREFIX}}_L001_R1_001.fastq*_Cut_0.fastq.gz {output.fq1}
        mv ${{OUTPREFIX}}_L001_RN_001.fastq*_MBC_0.txt.gz {output.mbc}
        mv ${{OUTPREFIX}}_L001_RN_001.fastq*_STATS_0.properties {output.properties}
        """

rule picard_MergeBamAlignment_dev_for_locatit_without_mbc:
    input:
        "out/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60.bam"

rule picard_FastqToSam_pe_dev_for_locatit_without_mbc:
    input:
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz"
    output:
        bam = "out/{tool}{extra}/{filler}.bam"
    wildcard_constraints:
        tool = "picard/FastqToSam"
    shell:
        """
        java -jar picard.jar FastqToSam \
            F1={input.fq1} \
            F2={input.fq2} \
            O={output.bam} \
            SM=sample001 \
            RG=rg0013
        """

rule agent_locatit_mbc:
    """
    Aim:
        This version of Locatit requires MBC file, and is for example required for RNA-seq in XTHS2 mode (Single-strand consensus)
    Doc:
        https://www.agilent.com/cs/library/software/public/AGeNT%20ReadMe.pdf
    Test:
        out/agent/locatit_mbc_-i_-R/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60.bam
    """
    input:
        #bam = "out/star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/agent/trim_-v2/ln/updir/mw/inp/fastq/2021_RNAseq_NECKER_spicuglia/fastq/MOLT4_S60.bam",
        # star/pe_fastq.gz_to_bam_staridx-GRCh38-ensembl_gtf-GRCh38-ensembl/
        # 
        bam = "out/{filler_align}{filler_trim}.bam",
        mbc = "out/{filler_trim}.txt.gz",
        agent = "out/agent/agent/agent.sh",
        locatit = "out/agent/agent/lib/locatit-2.0.5.jar"
    output:
        bam = "out/{tool}{extra}/{filler_align}{filler_trim}.bam"
    params:
        extra = params_extra,
        memory="128G" # default not working for test bam, this is the first other tested value
    wildcard_constraints:
        tool = "agent/locatit_mbc",
        filler_trim = "agent/trim.*"
    shell:
        """
        java -Xmx{params.memory} -jar {input.locatit} {params.extra} -o {output.bam} {input.bam} {input.mbc}

        # {input.agent} -Xmx{params.memory} locatit {params.extra} -o {output.bam} {input.bam} {input.mbc}
        """

rule get_agent:
    """
    https://explore.agilent.com/AGeNT-Software-Download-Form-TY
    """
    output:
        agent = "out/agent/agent/agent.sh",
        locatit = "out/agent/agent/lib/locatit-2.0.5.jar"
    shell:
        """
        OUTDIR=out/agent
        mkdir -p $OUTDIR
        cd $OUTDIR
        wget 'https://dt4ei3l3hxs7z.cloudfront.net/?elqTrackId=30b3c5b8c3bd44f7b3a01b66ab2a30a5&elqaid=3928&elqat=2' --output-document 'AGeNT_2.0.5.zip'
        unzip AGeNT_2.0.5.zip
        """
