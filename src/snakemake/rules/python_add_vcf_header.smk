rule python_add_vcf_header:
    """
    Created:
        2024-07-17 10:46:42
    Aim:
        Add a header to a VCF file. This was created to fix "vcf" file from the SMaSH repository: 
        https://github.com/rbundschuh/SMaSH
        The header is not needed for the SMaSH script but is needed for bedtools intersect to recognize the file as VCF and not bed-like.
    Test:
        out/python/add_vcf_header/wget/https/raw.githubusercontent.com/rbundschuh/SMaSH/master/snps_GRCh38.vcf
    """
    input:
        vcf = "out/{filler}.vcf",
        py = "../mw-lib/src/python/add_vcf_header.py"
    output:
        vcf = "out/python/add_vcf_header/{filler}.vcf"
    log:
        "out/python/add_vcf_header/{filler}.log"
    shell:
        "{input.py} {input.vcf} {output.vcf} &> {output}.log"