rule cuffmerge:
    """
    Created:
        2017-03-26 19:46:54
    Aim:
        Merge together several cufflinks assemblies, making it easier to produce an assembly GTF file suitable for use with cuffdiff.
    Note:
        # Denis code:
        #ls -1 output/cufflinks/*/transcripts.gtf > output/cuffmerge/assembly.txt
        #cuffmerge -o output/cuffmerge -g {params.gtf} --keep-tmp -s {params.fa} -p 5 output/cuffmerge/assembly.txt &> {output}.log
    Test:
        out/cuffmerge/GRCm38_H2AL2/merged.gtf
    """
    input:
        ref_sequence="out/cat/assembly_ensembl/{index}.fa",
        ref_gtf="out/ln/annotation_ensembl/{index}.gtf",
        cufflinks_gtf=expand(
            "out/cufflinks/g-{{index}}_lt-fr-firststrand/star/pe_{{index}}_outFilterMultimapNmax-{{outFilterMultimapNmax}}/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-{stage}-H2AL2-{condition}-Rep{replicate}/transcripts.gtf",
            stage=["P","R","C"],
            condition=["WT","KO"],
            replicate=["1","2","3"]
            )
    output:
        gtf="out/cuffmerge/{index}_H2AL2_outFilterMultimapNmax-{outFilterMultimapNmax}/merged.gtf",
        txt="out/cuffmerge/{index}_H2AL2_outFilterMultimapNmax-{outFilterMultimapNmax}/assemblies.txt"
        #log="out/cuffmerge/{index}_H2AL2/log.txt"
    params:
        outdir="out/cuffmerge/{index}_H2AL2_outFilterMultimapNmax-{outFilterMultimapNmax}"
    wildcard_constraints:
        outFilterMultimapNmax="[0-9]+"
    threads:
        4
    conda:
        "../envs/cufflinks.yaml"
    shell:
        """
        # Preparing input file for cuffmerge
        echo "{input.cufflinks_gtf}" | tr ' ' '\n' > {output.txt}
        echo "{output.txt}"

        cuffmerge \
            -o {params.outdir} \
            --ref-gtf {input.ref_gtf} \
            --ref-sequence {input.ref_sequence} \
            --num-threads {threads} \
            {output.txt}
        """

