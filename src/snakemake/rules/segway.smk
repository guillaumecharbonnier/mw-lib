rule segway_quickstart_traindir:
    """
    Created:
        2017-09-07 15:10:39
    """
    input:
        genomedata="out/wget/http/pmgenomics.ca/hoffmanlab/proj/segway/2011/test.genomedata"
    output:
        train_tab="out/segway/quickstart/traindir/train.tab"
    params:
        traindir="out/segway/quickstart/traindir"
    conda:
        "../envs/segway.yaml"
    shell:
        "segway --num-labels=4 train {input.genomedata} {params.traindir}"

rule segway_quickstart_identify:
    """
    Created:
        2017-09-07 15:52:54
    """
    input:
        genomedata="out/wget/http/pmgenomics.ca/hoffmanlab/proj/segway/2011/test.genomedata",
        train_tab="out/segway/quickstart/traindir/train.tab"
    output:
        bed="out/segway/quickstart/identifydir/segway.bed.gz",
        layered_bed="out/segway/quickstart/identifydir/segway.layered.bed.gz"
    params:
        traindir="out/segway/quickstart/traindir",
        identifydir="out/segway/quickstart/identifydir"
    conda:
        "../envs/segway.yaml"
    shell:
        "segway identify {input.genomedata} {params.traindir} {params.identifydir}"
