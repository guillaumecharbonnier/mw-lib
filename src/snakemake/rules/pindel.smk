rule pindel:
    """
    TODO: look at transindel
    https://github.com/tommyau/indelseek
    https://github.com/NCGG-MGC/IMSindel
    https://github.com/ratan-lab/indelMINER
    gatk4

    Doc:
        http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
        https://www.biostars.org/p/314926/
        -c 5:135,402,841-135,402,862

        GRCh38 coordinates where TAL1 insertion is located
        chr1:47,208,640-47,283,746

        Another option is IMSindel, that according to the paper has good performance when compared to others like Pindel.

    Test:
        out/pindel/_-c_1:47,208,640-47,283,746_fa-genome-GRCh38/ln/updir/mw-tall/src/pindel/jurkat/done
    """
    input:
        #bam="out/{filler}.bam", # TODO: write function that parse first column of Txt to get bam dependencies
        txt="out/{filler}.txt",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf=touch("out/{tool}{extra}_{fa_genome_id}/{filler}/done")
    log:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/log"
    benchmark:
        "out/{tool}{extra}_{fa_genome_id}/{filler}/benchmark.tsv"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="pindel/"
    threads:
        MAX_THREADS
    conda:
        "../envs/pindel.yaml"
    shell:
        """
        pindel --number_of_threads {threads} -i {input.txt} -f {input.fa} -o `dirname {output}` &> {log}
        """

rule pindel2vcf:
    """
    TODO: look at transindel
    Doc:
        http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
    """
    input:
        #bam="out/{filler}.bam",
        done="out/{filler}/done",
        fa= lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        vcf="out/{tool}{extra}_{fa_genome_id}/{filler}.vcf"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="pindel2vcf/"
    #threads:
    #    MAX_THREADS
    conda:
        "../envs/pindel.yaml"
    shell:
        """
        pindel2vcf --pindel_output_root `dirname {input.done}` -r {input.fa} -o `dirname {output}`
        """

    
#pindel2vcf --pindel_output_root output/Sample8 -r /ReferenceMaterial/1000Genomes/human_g1k_v37.fasta -R 1000GenomesPilot-NCBI36 -d 20101123 --min_coverage 18 --het_cutoff 0.2 --hom_cutoff 0.8 --vcf output/Sample8.vcf
#bgzip -f output/Sample8.vcf
#tabix -f -p vcf output/Sample8.vcf.gz
#bcftools view -Ov --exclude-uncalled --min-ac=1 output/Sample8.vcf.gz > output/Sample8.filt.vcf
