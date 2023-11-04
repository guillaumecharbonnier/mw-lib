rule rm_tmp:
    """
    Created:
        2017-03-13 15:35:31
    Aim:
        Remove a list of files and directories that are neither mandatory nor consuming to reproduce.
    """
    shell:
        """
        rm -rf \
            out/sort \
            out/unzip \
            out/gunzip \
            out/tar \
            out/fastx_toolkit \
            out/java \
            #out/samtools
            #out/cat \

        """

