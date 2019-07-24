rule find_delete_uncompressed_heavy_files:
    """
    Created:
        2019-03-22 13:25:41
    Aim:
        Slim down output directory.
    """
    shell:
        """
        find out/ \
            -name '*.sam' -o \
            -name '*.fastq' \
            -delete
        """
