rule python_convert_csfastq_to_solid:
    """
    Created:
        2018-01-25 11:08:57
    Aim:
        Convert xlsx files to tsv for easier view in terminal, and use as input for Salva workflow
    Test:
        out/python/convert_csfastq_to_solid/ln/abspath/gpfs/tagc/home/sadouni/reads/fastq/Run_93/Run_92_WF_14_11_24_RNAseq_Necker_L05_ISP_F3.qual
    """
    input:
        script="../mw-lib/src/python/convert_csfastq_to_solid.py",
        fastq="out/{filler}.fastq"
    output:
        csfasta=temp("out/python/convert_csfastq_to_solid/{filler}.csfasta"),
        qual=temp("out/python/convert_csfastq_to_solid/{filler}.qual"),
    #conda:
    #    "../envs/openpyxl.yaml"
    shell:
        "{input.script} {input.fastq} {output.csfasta} {output.qual}"

