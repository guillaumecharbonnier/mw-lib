rule ucsc_hgGcPercent:
    """
    Created:
        2017-07-28 10:54:09
    Aim:
        Testing this tool to produce a GCpercent track for mm10 to use with Hilbert Curves.
    Test:
        out/ucsc/hgGcPercent/2bit-mm10.wig
    """
    input:
        twobit = lambda wildcards: config['ids'][wildcards.twobit_id]
    output:
        wig="out/{tool}{extra}/{twobit_id}.wig"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="ucsc/hgGcPercent"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "hgGcPercent {params.extra} -file={output.wig} {input.bit}"
        #"hgGcPercent -wigOut -doGaps -file=stdout -win=200 mm10 {input.bit}"

rule ucsc_hgGcPercent_bigwig:
    """
    Created:
        2017-07-28 10:54:09
    Aim:
        Testing this tool to produce a GCpercent track for mm10 to use with Hilbert Curves.
    """
    input:
        chromInfo="out/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt",
        bit="out/rsync/ucsc/goldenPath/mm10/bigZips/mm10.2bit"
        #bit="out/wget/http/hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit"
    output:
        bw="out/ucsc/hgGcPercent/mm10.bw"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "hgGcPercent -wigOut -doGaps -file=stdout -win=200 mm10 {input.bit} |"
        "wigToBigWig -clip stdin {input.chromInfo} {output.bw}"
