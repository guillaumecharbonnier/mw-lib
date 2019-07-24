rule gem_indexer:
    """
    Created:
        2017-07-27 13:59:53
    """
    input:
        gem_indexer="opt/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-indexer",
        fa="out/cat/assembly_ensembl/GRCm38.fa"
    output:
        gem="out/gem/indexer/GRCm38.gem",
        log="out/gem/indexer/GRCm38.log"
    params:
        outprefix="out/gem/indexer/GRCm38"
    shell:
        """
        export PATH=$PATH:opt/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
        mkdir -p {params.outdir}

        gem-indexer \
            -i {input.fa} \
            -o {params.outprefix}
        """

rule gem_mappability:
    """
    Created:
        2017-07-27 14:19:15
    """
    input:
        gem_indexer="opt/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-mappability",
        idx="out/gem/indexer/GRCm38.gem"
    output:
        map="out/gem/mappability/GRCm38.mappability"
    params:
        outprefix="out/gem/mappability/GRCm38"
    threads:
        16
    shell:
        """
        export PATH=$PATH:opt/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin

        gem-mappability \
            -I {input.idx} \
            -l 60 \
            -T {threads} \
            -o {params.outprefix}
        """

rule gem2wig:
    """
    Created:
        2017-07-28 15:26:43
    Aim:
        Convert gem file to wig.
    Note:
        long arguments do not work.
    """
    input:
        gem2wig="opt/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin/gem-2-wig",
        map="out/gem/mappability/GRCm38.mappability",
        idx="out/gem/indexer/GRCm38.gem"
    output:
        wig="out/gem/gem-2-wig/GRCm38.wig"
    params:
        outprefix="out/gem/gem-2-wig/GRCm38"
    shell:
        """
        export PATH=$PATH:opt/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin

        gem-2-wig \
            -I {input.idx} \
            -i {input.map} \
            -o {params.outprefix}
        """
