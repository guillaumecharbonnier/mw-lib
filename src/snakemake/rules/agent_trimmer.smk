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
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz",
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
        FQ1=$OUTDIR/`basename {input.fq1}`
        FQ2=$OUTDIR/`basename {input.fq2}`

        # Input fastq symlink are dereferenced in the output directory
        # else agent will put output files into symlink targets.
        cp -rL {input.fq1} $FQ1
        cp -rL {input.fq2} $FQ2

        OUTPREFIX=out/{wildcards.tool}{wildcards.extra}/{wildcards.filler}

        # Output files from previous run should be removed
        # as the regex `mv` below should match only one file each.
        rm -f \
            ${{OUTPREFIX}}_1.fastq*_Cut_0.fastq.gz \
            ${{OUTPREFIX}}_2.fastq*_Cut_0.fastq.gz \
            ${{OUTPREFIX}}_N*_MBC_0.txt.gz \
            ${{OUTPREFIX}}_N*_STATS_0.properties

        {input.agent} trim -fq1 $FQ1 -fq2 $FQ2 {params.extra}

        rm -f $FQ1 $FQ2

        # Timestamped output files from agent are simplified
        mv ${{OUTPREFIX}}_2.fastq*_Cut_0.fastq.gz {output.fq2}
        mv ${{OUTPREFIX}}_1.fastq*_Cut_0.fastq.gz {output.fq1}
        mv ${{OUTPREFIX}}_N*_MBC_0.txt.gz {output.mbc}
        mv ${{OUTPREFIX}}_N*_STATS_0.properties {output.properties}
        """



rule agent_trim_2_lanes_pe:
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
        tool = "agent/trim_2_lanes"
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
