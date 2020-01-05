rule bbmap_readlength:
    """
    Created:
        2019-12-02 01:30:30
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/bbmap/readlength_csfasta.gz/ln/alias/experiments/solid_rnaseq_tall/ALL-T-2.10_12_F5.txt
    """
    input:
        "out/{filler}.{ext}",
    output:
        txt="out/{tool}_{ext}{extra}/{filler}.txt"
    log:
        "out/{tool}_{ext}{extra}/{filler}.log"
    benchmark:
        "out/{tool}_{ext}{extra}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "bbmap/readlength",
        ext = "fa|fq|sam|fastq.gz|csfasta.gz" #ADD HERE OTHER SUPPORTED EXTENSIONS I CURRENTLY HAVE NO INTEREST IN
    conda:
        "../envs/bbmap.yaml"
    shell:
        "readlength.sh in={input} out={output.txt} {params.extra} &> {log}"

