rule bismark_genome_preparation_fa:
    input:
        ref = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        directory("out/bismark/genome_preparation/{fa_genome_id}/Bisulfite_Genome"),
        ref = "out/bismark/genome_preparation/{fa_genome_id}/genome.fa"
        # "out/bismark/genome_preparation/{fa_genome_id}"
        # directory("out/bismark/genome_preparation/{fa_genome_id}/Bisulfite_Genome")
    log:
        "out/bismark/genome_preparation/{fa_genome_id}/log"
    params:
        ""  # optional params string
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        outdir=$(dirname {output.ref})
        echo $outdir
        ln -srf {input.ref} $outdir/genome.fa
        cd $outdir
        bismark_genome_preparation --verbose  .
        """
