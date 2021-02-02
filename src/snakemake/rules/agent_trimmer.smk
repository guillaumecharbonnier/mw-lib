rule agent_trim_pe:
    """
    Doc:
        https://www.agilent.com/cs/library/software/public/AGeNT%20ReadMe.pdf
        https://www.agilent.com/cs/library/software/public/AGeNTBestPractices.pdf
        
    Note:
        -out_loc `dirname {output.fq1}`
        does not work properly as output files are put inside input directory instead...
    Test:
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
        mbc = "out/{tool}{extra}/{filler}_txt.gz",
        properties = "out/{tool}{extra}/{filler}.properties",
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "agent/trim"
    shell:
        """
        rm -f \
            out/{wildcards.filler}_L001_R[1-2]_001.fastq[0-9]+_Cut_0.fastq.gz \
            out/{wildcards.filler}_L001_RN_001.fastq[0-9]+_MBC_0.txt.gz \
            out/{wildcards.filler}_L001_RN_001.fastq[0-9]+_STATS_0.properties

        {input.agent} trim -fq1 {input.fq1l1},{input.fq1l2} -fq2 {input.fq2l1},{input.fq2l2} {params.extra}

        mv out/{wildcards.filler}_L001_R1_001.fastq[0-9]+_Cut_0.fastq.gz {output.fq1}
        mv out/{wildcards.filler}_L001_R2_001.fastq[0-9]+_Cut_0.fastq.gz {output.fq2}
        mv out/{wildcards.filler}_L001_RN_001.fastq[0-9]+_MBC_0.txt.gz {output.mbc}
        mv out/{wildcards.filler}_L001_RN_001.fastq[0-9]+_STATS_0.properties {output.properties}
        """

rule get_agent:
    """
    https://explore.agilent.com/AGeNT-Software-Download-Form-TY
    """
    output:
        agent = "out/agent/agent/agent.sh"
    shell:
        """
        OUTDIR=out/agent
        mkdir -p $OUTDIR
        cd $OUTDIR
        wget 'https://dt4ei3l3hxs7z.cloudfront.net/?elqTrackId=30b3c5b8c3bd44f7b3a01b66ab2a30a5&elqaid=3928&elqat=2' --output-document 'AGeNT_2.0.5.zip'
        unzip AGeNT_2.0.5.zip
        """