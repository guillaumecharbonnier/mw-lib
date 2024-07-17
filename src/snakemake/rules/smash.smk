rule smash_dir:
    """
    Aim:
        First used to select subpopulations from bam.
        Trying this:
        https://www.biostars.org/p/225198/
    Note:
    Test:
        out/smash/vcf-SMaSH-snps-GRCh38/ln/alias/sst/by_run/run356/GRCh38/bam/pval_out.txt

    """
    input:
        # bam_list = lambda wildcards: config["bam_list"][wildcards.filler],
        # BAMs must end in .bam and be indexed
        vcf = lambda wildcards: eval(mwconf['ids'][wildcards.vcf_id]),
        py = "out/wget/https/raw.githubusercontent.com/rbundschuh/SMaSH/master/SMaSH.py"
    output:
        # TODO: Look for multiqc_dir rule to do the same
        tsv = "out/{tool}{extra}{vcf_id}/{filler}/pval_out.txt"
    params:
        indir = "out/{filler}",
        outdir = "out/{tool}{extra}{vcf_id}/{filler}"
    log:
        "out/{tool}{extra}{vcf_id}/{filler}.log"
    benchmark:
        "out/{tool}{extra}{vcf_id}/{filler}.benchmark.tsv"
    params:
        extra = params_extra
    threads:
        MAX_THREADS
    wildcard_constraints:
        tool = "smash/"
    conda:
        "../envs/smash.yaml"
    shell:
        """
        (
        mkdir -p {params.outdir}

        for ext in bam bam.bai
        do
            for file in "{params.indir}/"*.$ext
            do
                ln -srf "$file" "{params.outdir}"
            done
        done
        cd {params.outdir}

        python {WDIR}/{input.py} -i {WDIR}/{input.vcf} -bam ALL
        ) &> {log}
        """
