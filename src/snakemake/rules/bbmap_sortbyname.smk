rule bbmap_sortbyname:
    """
    Created:
        2017-11-12 22:22:20
    Aim:
        Try to sort genome.fa into an order that do not trigger error for yumei script.
    Test:
        out/bbmap/sortbyname/tar/xvzf_igenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
    """
    input:
        fa="out/{filler}.fa"
    output:
        fa="out/bbmap/sortbyname/{filler}.fa"
    conda:
        "../envs/bbmap.yaml"
    shell:
        "sortbyname.sh in={input.fa} out={output.fa}"

