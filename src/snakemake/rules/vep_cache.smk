rule vep_cache:
    """
    Derived from: https://gist.github.com/ckandoth/61c65ba96b011f286220fa4832ad2bc0
    """
    # output:
    #     "vep_cache/{build}.tar.gz"
    shell:
        """
        mkdir out/vep_cache
        cd out/vep_cache
        wget ftp://ftp.ebi.ac.uk/ensemblorg/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh37.tar.gz
        tar -zxf homo_sapiens_vep_102_GRCh37.tar.gz
        wget ftp://ftp.ebi.ac.uk/ensemblorg/pub/grch37/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
        
        gunzip Homo_sapiens.GRCh37.dna.toplevel.fa.gz
        bgzip -i Homo_sapiens.GRCh37.dna.toplevel.fa
        samtools faidx Homo_sapiens.GRCh37.dna.toplevel.fa.gz
        """
