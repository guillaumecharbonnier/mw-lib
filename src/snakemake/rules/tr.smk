rule tr_fasta_to_uppercase:
    """
    Created:
        2017-02-07 16:33:01
    Aim:
        In the assemblies there are sometime lowercase nucleotides. I have weird peak-motids results so maybe it is caused by the presence of the 'a','t','c','g' instead of 'A','T','C','G'.
    Note:
        Use of tr instead of awk because I suppose it is faster but not checked.
    Test:
       out/tr/fasta_to_uppercase/bedtools/getfasta/assemblymm9/bedtools/intersect_extract_ss_in_nuc/test2_nuc.fa 
       out/tr/fasta_to_uppercase/bedtools/getfasta/assemblymm9/bedtools/intersect_extract_ss_in_nuc/test2_ss.fa 
    """
    input:
        fa="out/{filler}.fa"
    output:
        fa="out/tr/fasta_to_uppercase/{filler}.fa"
    shell:
        """
        cat {input.fa} | tr '[atcg]' '[ATCG]' > {output.fa}
        """

rule tr_bed_to_gprofiler:
    """
    Created:
        2018-11-25 19:44:12
    Aim:
    Test:
        out/tr/bed_to_gprofiler/ln/alias/experiments/encode_broad_h3k4me3/Cancer_K562.txt
    """
    input:
        bed="out/{filler}.bed"
    output:
        txt="out/tr/bed_to_gprofiler/{filler}.txt"
    shell:
        """
        cut -f1-3 {input.bed} | tr '\\t' ':' > {output.txt}
        """
