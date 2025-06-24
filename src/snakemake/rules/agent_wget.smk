rule agent_wget:
    """
    Aim:
        Agent is currently only available on Agilent website
        https://explore.agilent.com/AGeNT-Software-Download-Form-TY
    """
    output:
        agent = "out/agent/agent/agent.sh",
        # Note locatit was removed in later versions of AGeNT. Currently commented but it would require adjustment to downstream RNA-Seq analyses.
        # locatit = "out/agent/agent/lib/locatit-2.0.5.jar"
    shell:
        """
        OUTDIR=out/agent
        mkdir -p $OUTDIR
        cd $OUTDIR
        # TODO: version is hardcoded in name but not the target url.
        # Need to find a specific URL to get the needed version.
        # This is actually not the 2.0.5 version.
        # The downstream rules expect these libs:
            # -rw-rw-r-- 1 gcharbonnier shinyusers  20K Aug 26  2022 Agent-3.0.6.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers 4,5M Mar  9  2022 creak-1.0.5.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers 119K Nov 18  2021 junit-3.8.1.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers 6,3M Sep 25  2020 locatit-2.0.5.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers  22M Aug 26  2022 trimmer-3.0.5.jar
        # But the current download leads to these libs:
            # -rw-rw-r-- 1 gcharbonnier shinyusers  16K Jun 28  2024 Agent-3.1.2.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers 4,3M Apr 24  2024 creak-1.1.1.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers 119K Nov 18  2021 junit-3.8.1.jar
            # -rw-rw-r-- 1 gcharbonnier shinyusers  22M Jun 28  2024 trimmer-3.1.2.jar
        wget 'https://dt4ei3l3hxs7z.cloudfront.net/?elqTrackId=30b3c5b8c3bd44f7b3a01b66ab2a30a5&elqaid=3928&elqat=2' --output-document 'AGeNT_2.0.5.zip'
        unzip AGeNT_2.0.5.zip
        chmod +x agent/agent.sh
        """
