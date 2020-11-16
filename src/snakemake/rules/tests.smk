"""
Extract tests string
sed -n '/Test:/,/""/p' *.rules | sed -n '/Test:/!p' | sed -n '/""/!p'
"""

rule tests_light:
    input:
        "out/wget/ftp/ftp.ensembl.org/robots.txt"


rule tests:
    input:
        "out/ucsc/hgGcPercent/2bit-mm10.wig"
