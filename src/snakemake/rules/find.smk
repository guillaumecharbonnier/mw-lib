rule find_delete_uncompressed_heavy_files:
    """
    Created:
        2019-03-22 13:25:41
    Aim:
        Slim down output directory.
    """
    log: "out/find/delete_uncompressed_heavy_files.log"
    shell:
        """
        find out/ \
            -name '*.sam' -o \
            -name '*.fastq' \
            -delete &> {log}
        """

rule find_delete_broken_symlinks:
    """
    Created:
        2019-03-22 13:25:41
    Aim:
        After using a rule that remove some output files, It may be needed to remove broken symlinks to avoid Snakemake errors.
        This rule create a txt file with the list of broken symlinks and then remove them.
        The txt file can be used to check the deleted files.

        WIP: This rule does not work as expected.
        I actually found a way to remove the related snakemake error by specifyin --default-resources 'mem_mb=1000' 'disk_mb=1000'
        Doing so avoid the default call to "input.size_mb" that throw the error for broken symlinks.
    """
    log:
        "out/find/delete_broken_symlinks.log"
    shell:
        """
            find out/ -xtype l -print -delete &> {log}
        """
