rule pylint:
    """
    Created:
        2019-03-04 16:07:57
    Aim:
    Note:
        Currently not working
    Test:
        out/pylint/ln/updir/mw-sst/src/snakemake/functions/mapper2.snake.txt"
    """
    input:
        "out/{filler}"
    log:
        "out/pylint/{filler}.txt"
    conda:
        "../envs/pylint.yaml"
    shell:
        "pylint {input} > {output}"

rule pyreverse:
    """
    Test:
        out/pyreverse/ln/updir/mw-sst/src/snakemake/functions/mapper.dot
    """
    input:
        py  = "out/{filler}.py"
    output:
        py  = "out/pyreverse/{filler}.py",
        png = "out/pyreverse/{filler}.png"
        #dot = "out/pyreverse/{filler}.dot"
    conda:
        "../envs/pylint.yaml"
    shell:
        "ln -srf {input.py} {output.py}; "
        "cd `dirname {output.py}`; "
        "pyreverse -o png `basename {output.py}`"
