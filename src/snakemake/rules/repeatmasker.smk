rule repeatmasker_test1:
    """
    out/repeatmasker/test2_-species_mouse_-a_-html/rnaspades/pe_test_h2al2/transcripts/transcripts.fasta.out

    $ ll out/repeatmasker/test1_-species_mouse_-a_-html/rnaspades/pe_test_h2al2/transcripts/
    -rw-rw-r-- 1 Charbonnier thymus  82M 31 janv. 01:10 transcripts.fasta.cat.gz
    -rw-rw-r-- 1 Charbonnier thymus 558M 31 janv. 01:10 transcripts.fasta.out.html
    -rw-rw-r-- 1 Charbonnier thymus  223 31 janv. 01:10 transcripts.fasta.alert
    -rw-rw-r-- 1 Charbonnier thymus 340M 31 janv. 01:10 transcripts.fasta.align
    -rw-rw-r-- 1 Charbonnier thymus  49M 31 janv. 01:10 transcripts.fasta.out
    -rw-rw-r-- 1 Charbonnier thymus 201M 31 janv. 01:10 transcripts.fasta.masked
    -rw-rw-r-- 1 Charbonnier thymus 2,0K 31 janv. 01:10 transcripts.fasta.tbl


    """
    input:
        fasta="out/{filler}/{filename}.fasta"
    output:
        expand("out/{{tool}}{{extra}}/{{filler}}/{{filename}}/{{filename}}.fasta.{ext}", ext=["out"])
    params:
        extra = params_extra
    wildcard_constraints:
        tool="repeatmasker/test2",
        filename="[a-zA-Z0-9-_]+"
    conda:
        "../envs/repeatmasker.yaml"
    threads:
        16
    shell:
        """
        RepeatMasker -parallel {threads} {params.extra} --dir `dirname {output}` {input.fasta}
        """
