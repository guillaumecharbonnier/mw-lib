rule multiqc_dir:
    """
    Created:
        2017-03-23 11:20:51
    Aim:
        Aggregate results from bioinformatics analyses across many samples into a single report
        This version crawls available log file in the *filler* directory.
    Note:
        Remove first output directory because else output html has a different name (with suffix)
    Test:
        out/multiqc/dir/out/multiqc_report.html
        out/multiqc/dir/ln/alias/sst/by_customer/Salva_Vahid/logs/multiqc_report.html
    """
    output:
        "out/{tool}{extra}/{filler}multiqc_report.html"
    log:
        "out/{tool}{extra}/{filler}multiqc_report.log"
    benchmark:
        "out/{tool}{extra}/{filler}multiqc_report.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "multiqc/dir"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "rm -rf `dirname {output}`/multiqc_data; "
        "multiqc --outdir `dirname {output}` {params.extra} out/{wildcards.filler} &> {log}"

rule create_multiqc_filelist:
    """
    Test:
        out/create_multiqc_filelist/multiqc-sst.txt
    """
    input:
        #dependencies = lambda wildcards: eval(str(config['ids'][wildcards.multiqc_id]))
    output:
        filelist = "out/create_multiqc_filelist/{multiqc_id}.txt"
    run:
        with open(output.filelist, "w") as f:
            #f.writelines(map(lambda s: s + '\n', lines))
            for item in config['ids'][wildcards.multiqc_id]:
                f.write("%s\n" % item)

rule multiqc_with_requirements:
    """
    Created:
        2019-03-06 10:44:18
    Aim:
        Aggregate results from bioinformatics analyses across many samples into a single report.
        This version requires the log files to be used as input.
    Test:
        out/multiqc/req_multiqc-test1/multiqc_report.html
        out/multiqc/req_multiqc-sst/multiqc_report.html
    """
    input:
        #dependencies = lambda wildcards: eval(config['ids'][wildcards.multiqc_id])
        "out/create_multiqc_filelist/{multiqc_id}.txt"
    output:
        report   = "out/{tool}{extra}_{multiqc_id}/multiqc_report.html"
    log:
        "out/{tool}{extra}_{multiqc_id}/multiqc_report.log"
    benchmark:
        "out/{tool}{extra}_{multiqc_id}/multiqc_report.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool = "multiqc/req"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --outdir `dirname {output}` {params.extra} --file-list {input}"
        #"(multiqc --force --outdir `dirname {output}` {params.extra} {input}) &> {log}"

