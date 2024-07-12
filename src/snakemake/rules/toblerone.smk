
# Load the data
toblerone_test_data = pandas.read_csv('../mw-lib/src/snakemake/rules/toblerone.tsv', sep='\t')

# Get the 'Patient' column as a list
toblerone_test_cdna = toblerone_test_data['Patient'].tolist()

rule toblerone_tests_LEF1:
    input:
        expand(
            "out/toblerone/map_fa-genome-hg38_LEF1/agent/trim_-v2/ln/alias/sst/all_samples/fastq/{cdna}_RNA_XT_HS2.csv",
            cdna = toblerone_test_cdna
        )


rule toblerone_map_extra:
    """
    Created:
        2024-03-29 20:01:56
    Test:
        out/toblerone/map_fa-genome-hg38_IKZF1/agent/trim_-v2/ln/alias/sst/all_samples/fastq/113281_PICCL.csv
        out/toblerone/map_fa-genome-hg38_IKZF1/agent/trim_-v2/ln/alias/sst/all_samples/fastq/887_RNA_XT_HS2.csv
        out/toblerone/map_fa-genome-hg38_LEF1/agent/trim_-v2/ln/alias/sst/all_samples/fastq/887_RNA_XT_HS2.csv


    """
    input:
        fq1 = "out/{filler}_1.fastq.gz",
        fq2 = "out/{filler}_2.fastq.gz",
        tinyt = "out/wget/https/github.com/Oshlack/Toblerone/releases/download/v0.0.9/tinyt_amd64",
        tidx = "out/toblerone/build_index_{fa_genome_id}_{gene_name}/toblerone_transcriptome.tidx"
    output:
        csv = "out/{tool}_{fa_genome_id}_{gene_name}/{filler}.csv"
    wildcard_constraints:
        tool="toblerone/map",
        gene_name="[A-Za-z0-9]+",
        fa_genome_id = "fa-genome-hg38",
    conda:
        "../envs/toblerone.yaml"
    threads:
        8
    shell:
        """
        # {input.tinyt} map -o {output.csv} -n {threads} -i {input.tidx} {input.fq1} {input.fq2}
        {input.tinyt} map -o {output.csv} -n {threads} -i {input.tidx} <(gunzip -c {input.fq1}) <(gunzip -c {input.fq2})
        """

rule toblerone_build_index:
    """
    https://github.com/Oshlack/Toblerone

    Test:
        out/toblerone/build_index_fa-genome-hg38_TP73/ucsc.bed12
        out/toblerone/build_index_fa-genome-hg38_IKZF1/ucsc.bed12
    """
    input:
        # gtf   = lambda wildcards: eval(mwconf['ids'][wildcards.gtf_id])
        retrievebedscript = "../mw-lib/src/python/retrieve_bed12.py",
        # script = "out/wget/https/raw.githubusercontent.com/Oshlack/Toblerone/master/scripts/create_bedfiles.py",
        script = "../mw-lib/src/python/create_bed_of_all_exon_deletions.py",
        tinyt = "out/wget/https/github.com/Oshlack/Toblerone/releases/download/v0.0.9/tinyt_amd64",
        fa_genome = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id]),
    output:
        ucsc_bed12 = "out/toblerone/build_index_{fa_genome_id}_{gene_name}/ucsc.bed12",
        combined_bed12 = "out/toblerone/build_index_{fa_genome_id}_{gene_name}/combined.bed12",
        toblerone_fa = "out/toblerone/build_index_{fa_genome_id}_{gene_name}/toblerone_transcriptome.fa",
        tidx = "out/toblerone/build_index_{fa_genome_id}_{gene_name}/toblerone_transcriptome.tidx"
    params:
        tmpdir = "out/toblerone/build_index_{fa_genome_id}_{gene_name}/tmp",
        # genome should be fa_genome_id without fa-genome- prefix
        genome_id = lambda wildcards: wildcards.fa_genome_id.split("fa-genome-")[1]
    conda:
        # "../envs/toblerone.yaml"
        "../envs/toblerone.yaml"
    wildcard_constraints:
        # Note that only genome id recognized by UCSC can be used here, and not all those defined in mwconf['ids']
        fa_genome_id = "fa-genome-hg38",
    shell:
        """
        mkdir -p {params.tmpdir}

        python {input.retrievebedscript} {wildcards.gene_name} {params.genome_id} {output.ucsc_bed12}

        python {input.script} {output.ucsc_bed12} {params.tmpdir}
        find {params.tmpdir}  -name "*.bed" -type f -exec cat {{}} + > {output.combined_bed12}
        rm -rf {params.tmpdir}

        bedtools getfasta -fi {input.fa_genome} -bed {output.combined_bed12} -split -name | \
            sed 's/::.*//' | \
            awk -v gene_name="{wildcards.gene_name}" '/^>/ {{print $0" gene="gene_name; next}} {{print}}' > {output.toblerone_fa}

        chmod +x {input.tinyt}
        {input.tinyt} index -i {output.tidx} {output.toblerone_fa}
        """


