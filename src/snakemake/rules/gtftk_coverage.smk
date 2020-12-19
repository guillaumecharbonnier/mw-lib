rule gtftk_coverage_extra_chrominfo_bed_bw:
    """
    Created:
        2017-03-07 15:47:19
    Aim:
        Takes a GTF as input to compute bigwig coverage in regions of interest (promoter, transcript body, intron, intron_by_tx, tts…) or a BED6 to focus on user-defined regions.
    Note:
        Pseudocount forced to 0 because of a strange bug with default pseudocount=1:
            https://github.com/dputhier/gtftk/issues/125
    Doc:
        https://dputhier.github.io/pygtftk/coverage.html#coverage
    Test:
        out/gtftk/coverage_u-100000_d-100000_p-0_w-1000_n-10_x_chrominfo-hg38-main-chr_bed-hg38-bs-hypometh-thymus-union-all_bw-hg38-RNA-thymus.txt
        out/gtftk/coverage_-x_chrominfo-mm10_bed-mm10-ikE120-peak_bw-mm10-immgen-atac.txt
    """
    input:
        inputfile = lambda wildcards: eval(config['ids'][wildcards.bed_list_id]),
        chrominfo = lambda wildcards: eval(config['ids'][wildcards.chrominfo_id]),
        bw        = lambda wildcards: eval(config['ids'][wildcards.bw_list_id])
    output:
        txt="out/{tool}{extra}_{chrominfo_id}_{bed_list_id}_{bw_list_id}.txt"
    log:
        "out/{tool}{extra}_{chrominfo_id}_{bed_list_id}_{bw_list_id}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="gtftk/coverage"
    threads:
        MAX_THREADS
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk coverage {params} --inputfile {input.inputfile} --outputfile {output.txt} --chrom-info {input.chrominfo} --nb-proc {threads} {input.bw} &> {log}
        """

rule gtftk_coverage_extra_chrominfo_bed_single_bw:
    """
    Aim:
        gtftk coverage version with only one bw file as input.
    Doc:
        https://dputhier.github.io/pygtftk/coverage.html#coverage
    Test:
        out/gtftk/coverage_-x_chrominfo-hg19_bed-hg19-ec-sharp/ln/alias/sst/all_samples/hg19/bw/EC_H3K4me3.bw
    """
    input:
        inputfile = lambda wildcards: eval(config['ids'][wildcards.bed_list_id]),
        chrominfo = lambda wildcards: eval(config['ids'][wildcards.chrominfo_id]),
        bw        = out/{filler}.bw
    output:
        txt="out/{tool}{extra}_{chrominfo_id}_{bed_list_id}/{filler}.txt"
    log:
        "out/{tool}{extra}_{chrominfo_id}_{bed_list_id}/{filler}.log"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="gtftk/coverage"
    threads:
        MAX_THREADS
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        gtftk coverage {params} --inputfile {input.inputfile} --outputfile {output.txt} --chrom-info {input.chrominfo} --nb-proc {threads} {input.bw} &> {log}
        """



rule gtftk_coverage_extra_chrominfo_gtf_bw:
    """
    Created:
        2017-03-07 15:47:19
    Aim:
        Same as gtftk_coverage_extra_chrominfo_bed_bw but using a gtf as input.
    Test:
        out/gtftk/coverage_p-0_x_f-exon_chrominfo-hg38-main-chr_gtf-hg38-ensembl_bw-hg38-RNA-thymus.txt
    """
    input:
        gtf       = lambda wildcards: eval(config['ids'][wildcards.gtf_id]),
        chrominfo = lambda wildcards: eval(config['ids'][wildcards.chrominfo_id]),
        bw        = lambda wildcards: eval(config['ids'][wildcards.bw_list_id])
    output:
        txt="out/{tool}{extra}_chrominfo-{chrominfo_id}_gtf-{gtf_id}_{bw_list_id}.txt"
    log:
        "out/{tool}{extra}_chrominfo-{chrominfo_id}_gtf-{gtf_id}_{bw_list_id}.log"
    threads:
        16
    wildcard_constraints:
        tool="gtftk/coverage"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        """
        (
        EXTRA=`echo {wildcards.extra} | sed -e 's/-/ /g' -e 's/_/ -/g'`
        gtftk coverage $EXTRA --inputfile {input.gtf} --outputfile {output.txt} --chrom-info {input.chrominfo} --nb-proc {threads} {input.bw}
        ) &> {log}
        """
