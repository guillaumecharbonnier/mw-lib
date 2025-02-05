
rule bsgenova_test:
    """
    Test:
        out/bsgenova/fa-genome-hg19/samtools/index/samtools/sort/fumi_tools/dedup_--paired/samtools/sort/bismark/pe_--pbat_fa-genome-hg19/trim_galore/pe_--non_directional_--rrbs_--length_15_--fastqc/ln/updir/mw-tall-rrbs/inp/bcl/240521_A00680_0489_BHVC2MDMXY/data-isilon/raw-data/novaseq_imagine/240521_A00680_0489_BHVC2MDMXY/Data/Intensities/BaseCalls/fumi_tools_demultiplex/638_GEOVAS.vcf.gz
    """
    input:
        bsextractor = "out/wget/https/raw.githubusercontent.com/hippo-yf/bsgenova/refs/heads/main/bsextractor.py",
        bsgenova = "out/wget/https/raw.githubusercontent.com/hippo-yf/bsgenova/refs/heads/main/bsgenova.py",
        bam = "out/{filler}.bam",
        # bsextractor expect the bai,
        bai = "out/{filler}.bam.bai",
        fa = lambda wildcards: eval(mwconf['ids'][wildcards.fa_genome_id])
    output:
        vcf = "out/bsgenova/{fa_genome_id}/{filler}.vcf.gz",
        snv = "out/bsgenova/{fa_genome_id}/{filler}.snv.gz",
        bed = "out/bsgenova/{fa_genome_id}/{filler}.bed.gz"
    params:
        prefix = "out/bsgenova/{fa_genome_id}/{filler}"
    conda:
        "../envs/bsgenova.yaml"
    shell:
        """
        python {input.bsextractor} --base-quality 20 --read-quality 20 --bam-file {input.bam} --reference-genome {input.fa} --output-bed {output.bed} --output-atcgmap - | python {input.bsgenova} -o {params.prefix}
        """