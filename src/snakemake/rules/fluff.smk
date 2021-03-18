rule fluff_heatmap:
    """
    Created:
        2019-02-18 14:20:07
    Doc:
        https://fluff.readthedocs.io/en/latest/commands.html#fluff-heatmap
    Test:
        out/fluff/heatmap_-C_k_-k_5/bed-mm10-test-srr-peaks_bam-mm10-test-srr.png
        out/fluff/heatmap_-C_k_-k_5_-g_-M_Pearson/bed-mm10-test-srr-peaks_bam-mm10-test-srr.png
    """
    input:
        bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id]),
        d   = lambda wildcards: eval(mwconf['ids'][wildcards.d_id])
    output:
        png = "out/{tool}{extra}/{bed_list_id}_{d_id}.{ext}",
        bed = "out/{tool}{extra}/{bed_list_id}_{d_id}.{ext}_clusters.bed",
        cnt = "out/{tool}{extra}/{bed_list_id}_{d_id}.{ext}_readCounts.txt"
    log:      "out/{tool}{extra}/{bed_list_id}_{d_id}.{ext}.log"
    benchmark:"out/{tool}{extra}/{bed_list_id}_{d_id}.{ext}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "fluff/heatmap",
        ext  = "png|pdf"
    conda:
        "../envs/fluff.yaml"
    threads:
        MAX_THREADS
    shell:
        "fluff heatmap {params} -f {input.bed} -d {input.d} -o {output.png} -P {threads} &> {log}"

rule fluff_bandplot:
    """
    Created:
        2019-02-18 14:20:07
    Doc:
        https://fluff.readthedocs.io/en/latest/commands.html#fluff-bandplot
        TODO: count should be used by default instead of d.
    Test:
        #out/fluff/bandplot/bed-mm10-test-srr-peaks-fluff-clusts_bam-mm10-test-srr.png
        out/fluff/bandplot_bam-mm10-test-srr/fluff/heatmap_-C_k_-k_5/bed-mm10-test-srr-peaks_bam-mm10-test-srr.png_clusters.png
    """
    input:
        #bed = lambda wildcards: eval(mwconf['ids'][wildcards.bed_list_id]),
        bed = "out/{filler}.bed",
        d   = lambda wildcards: eval(mwconf['ids'][wildcards.d_id])
    output:
        #png="out/{tool}{extra}/{bed_list_id}_{d_id}.png"
        img = "out/{tool}{extra}_{d_id}/{filler}.{ext}"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "fluff/bandplot",
        d_id = "[a-zA-Z0-9-_]+",
        ext = "png|pdf"
    conda:
        "../envs/fluff.yaml"
    threads:
        1
    shell:
        "fluff bandplot {params} -f {input.bed} -d {input.d} -o {output.img}"

rule fluff_profile:
    """
    Created:
        2019-02-18 14:20:07
    Doc:
        https://fluff.readthedocs.io/en/latest/commands.html#fluff-profile
    Test:
        out/fluff/profile_-i_chr2:112244957-112260300/bam-mm10-test-srr.png
    """
    input:
        d = lambda wildcards: eval(mwconf['ids'][wildcards.d_id])
    output:
        png="out/{tool}{extra}/{d_id}.png"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "fluff/profile"
    conda:
        "../envs/fluff.yaml"
    threads:
        1
    shell:
        "fluff profile {params} -d {input.d} -o {output.png}"

