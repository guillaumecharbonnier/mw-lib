rule agent_wget:
    """
    Aim:
        Agent is currently only available on Agilent website
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
