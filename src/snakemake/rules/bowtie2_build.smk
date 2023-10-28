rule bowtie2_build_tests:
    input:
        "out/bowtie2-build/fa-genome-GRCh38-r95-chr22.1.bt2",
        "out/bowtie2-build/fa-genome-GRCh38-r95-chr22.1.bt2l"

rule bowtie2_build:
    """
    Aim:
        Adjust bowtie2 build rule to allow the use of the cache directive.
    Doc:
        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    """
    input:
        ref = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        multiext(
            "out/bowtie2-build/{fa_genome_id}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "out/bowtie2-build/{fa_genome_id}.log"
    cache:
        "all"
    params:
        fa_genome_id = "{fa_genome_id}"
    threads:
        MAX_THREADS
    wrapper:
        "v1.28.0/bio/bowtie2/build"

rule bowtie2_build_large:
    input:
        ref = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id]),
    output:
        multiext(
            "out/bowtie2-build/{fa_genome_id}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    log:
        "out/bowtie2-build/{fa_genome_id}_large.log"
    cache:
        "all"
    params:
        fa_genome_id = "{fa_genome_id}",
        extra="--large-index",  # optional parameters
    threads:
        MAX_THREADS
    wrapper:
        "v1.28.0/bio/bowtie2/build"


