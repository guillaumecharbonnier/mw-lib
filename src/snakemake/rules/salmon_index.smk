rule salmon_index:
    """
    Test:
        out/salmon/index/wget/http/ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all
    """
    input:
        "out/{filler}.fa.gz"
    output:
        directory("out/salmon/index/{filler}")
    params:
        filler = "{filler}"
    cache:
        "omit-software"
    conda:
        "../envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i {output}"


rule salmon_index_refseq_with_decoys:
    """
    """
    input:
        rna = expand(
            "out/wget/https/ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.{number}.rna.fna.gz",
            number=range(1, 17)  # Specify the range of numbers from 1 to 16
        )
        dna = "out/wget/ftp/ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    output:
        directory("out/salmon/index_refseq_with_decoys/homo_sapiens")
    shell:
        """
        mkdir -p {output}
        cat {input} > {output}/transcriptome_genome.fa.gz

        zcat {input.dna} | grep "^>" | cut -d " " -f 1 | sed "s/>//g" > {output}/decoys.txt

        salmon index -t {output}/transcriptome_genome.fa.gz -d {output}/decoys.txt -i {output}

        rm -f {output}/transcriptome_genome.fa.gz
        """



rule salmon_index_ensembl_with_decoys:
    """
    Test:
        out/salmon/index_ensembl_with_decoys/release-102/fasta/homo_sapiens/Homo_sapiens.GRCh38
    """
    input:
        cdna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/cdna/{specie_assembly}.cdna.all.fa.gz",
        ncrna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/ncrna/{specie_assembly}.ncrna.fa.gz",
        dna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/dna/{specie_assembly}.dna.primary_assembly.fa.gz"
    output:
        directory("out/salmon/index_ensembl_with_decoys/{release_fasta_specie}/{specie_assembly}")
    params:
        release_fasta_specie = "{release_fasta_specie}",
        specie_assembly = "{specie_assembly}"
    cache:
        "omit-software"
    wildcard_constraints:
        release_fasta_specie="release-[0-9]+/fasta/[a-z_]+",
        specie_assembly="[A-Za-z0-9_.]+"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        mkdir -p {output}
        zcat {input} > {output}/transcriptome_genome.fa

        zcat {input.dna} | grep "^>" | cut -d " " -f 1 | sed "s/>//g" > {output}/decoys.txt

        salmon index -t {output}/transcriptome_genome.fa -d {output}/decoys.txt -i {output}

        rm -f {output}/transcriptome_genome.fa.gz
        """

rule salmon_index_ensembl_with_decoys_toplevel:
    """
    Aim:
        According to
        https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
        we should use the primary_assembly fasta for genome reference.
        It is however missing for some species, and equivalent to the available toplevel one.
    Test:
        out/salmon/index_ensembl_with_decoys_toplevel/release-102/fasta/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0
    """
    input:
        cdna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/cdna/{specie_assembly}.cdna.all.fa.gz",
        ncrna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/ncrna/{specie_assembly}.ncrna.fa.gz",
        dna = "out/wget/ftp/ftp.ensembl.org/pub/{release_fasta_specie}/dna/{specie_assembly}.dna.toplevel.fa.gz"
    output:
        directory("out/salmon/index_ensembl_with_decoys_toplevel/{release_fasta_specie}/{specie_assembly}")
    params:
        release_fasta_specie = "{release_fasta_specie}",
        specie_assembly = "{specie_assembly}"
    cache:
        "omit-software"
    wildcard_constraints:
        release_fasta_specie="release-[0-9]+/fasta/[a-z_]+",
        specie_assembly="[A-Za-z0-9_.]+"
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        mkdir -p {output}
        zcat {input} > {output}/transcriptome_genome.fa

        zcat {input.dna} | grep "^>" | cut -d " " -f 1 | sed "s/>//g" > {output}/decoys.txt

        salmon index -t {output}/transcriptome_genome.fa -d {output}/decoys.txt -i {output}

        rm -f {output}/transcriptome_genome.fa.gz
        """
