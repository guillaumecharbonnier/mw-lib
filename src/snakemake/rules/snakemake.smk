import os
from snakemake.utils import read_job_properties
"""
Explanation from:
https://bitbucket.org/snakemake/snakemake/wiki/Documentation

When executing a workflow on a cluster using the --cluster parameter (see below), Snakemake creates a job script for each job to execute. This script is then invoked using the provided cluster submission command (e.g. qsub). Sometimes you want to provide a custom wrapper for the cluster submission command that decides about additional parameters. As this might be based on properties of the job, Snakemake stores the job properties (e.g. rule name, threads, input files, params etc.) as JSON inside the job script. For convenience, there exists a parser function snakemake.utils.read_job_properties that can be used to access the properties.
"""

rule snakemake:
    """
    Created:
        2019-02-08 16:14:51
    Test:
        out/snakemake/stdout_--rulegraph_test.dot
        out/snakemake/stdout_--dag_broad-analysis.dot
    """
    output:
        txt = "out/{tool}{extra}"
    log:
              "out/{tool}{extra}.log"
    benchmark:
              "out/{tool}{extra}.benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="snakemake/stdout"
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake {params.extra} > {output.txt} 2> {log}"
