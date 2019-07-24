"""
Extract tests string
sed -n '/Test:/,/""/p' *.rules | sed -n '/Test:/!p' | sed -n '/""/!p'
"""

rule tests_light:
    input:

rule tests:
    input:
        "out/ucsc/hgGcPercent/2bit-mm10.wig"
