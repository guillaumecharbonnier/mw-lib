rule r_sort_computeMatrix_max_5_to_3_end:
    """
    Test 2016-12-13 13h17:
        input:
            "out/deepTools/computeMatrix/referencePoint/center/bamCoverage/mm9/MNase_2016_12_13/only_round_low_fuzz_for_article/200bp/all.txt.gz"
        output:
            "out/r/sortMatrixByMax5to3end/referencePoint/center/bamCoverage/mm9/MNase_2016_12_13/only_round_low_fuzz_for_article/200bp/all.txt.gz"
    """
    input:
        matrix="out/deepTools/computeMatrix/{id}.txt.gz",
        rscript="opt/miniconda/envs/r/bin/Rscript",
        code="code/r/script/sort_compute_matrix.R"
    output:
        matrix="out/r/sortMatrixByMax5to3end/{id}.txt.gz",
        gunzip=temp("out/r/sortMatrixByMax5to3end/{id}_gunzip_tmp.txt"),
        r_input=temp("out/r/sortMatrixByMax5to3end/{id}_input_tmp.txt"),
        r_output=temp("out/r/sortMatrixByMax5to3end/{id}_output_tmp.txt"),
        meta=temp("out/r/sortMatrixByMax5to3end/{id}_meta_tmp.txt"),
        data=temp("out/r/sortMatrixByMax5to3end/{id}_data_tmp.txt")
    shell:"""
    gunzip --stdout {input.matrix} > {output.gunzip}
    head -1 {output.gunzip} > {output.meta}
    tail -n+2 {output.gunzip} | sed 's/nan/0.000000/g' > {output.data}
    cat {output.meta} {output.data} > {output.r_input}
    
    {input.rscript} {input.code} -i {output.r_input} -o {output.r_output}
    cat {output.r_output} | gzip --stdout > {output.matrix}
    """
