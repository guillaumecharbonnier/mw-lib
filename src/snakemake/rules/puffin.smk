rule puffin:
    """
    Created: 2016-03-17 17h21
    Article:
    http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-S9-S11
    Documentation:
    http://alumni.cs.ucr.edu/~polishka/indexPuffin.html
    
    Puffin need bam for only one chromosome.
    Puffin need python2.7+.
    Usage
    
    Given that input reads are in input.bam, that contain only reads for particular chromosome
    ./bam2bedpe.sh input.bam > input.bed
    python Run.py input.bed
    The output will be printed in input.bed.nucs using next column format:
    <Position of the nucleosome center> <width of the peak> <confidence score> <"Fuzziness"> <Level of the curve that was used to detect nucleosome >
    
    Usage:
        "out/puffin/{filler}.bed.nucs"

        !!! WARNING !!! I noticed afterward that the bam2bedpe script does not output tab-delimited files... So it is not bed file. Check if Puffin really wants space-delimited files...
    """
    input:
        bam="out/samtools/view/split_bam_by_chr/{filler}.bam",
        bai="out/samtools/view/split_bam_by_chr/{filler}.bam.bai",
        python="opt/miniconda/envs/py27/bin/python",
        puffin="opt/puffin/Run.py",
        samtools="opt/samtools-1.3.1/samtools",
        bam2bedpe="opt/puffin/bam2bedpe.sh"
    output:
        bed="out/puffin/{filler}.bed",
        nucs="out/puffin/{filler}.bed.nucs"
    resources:
        ram=70 # Puffin seems to require huge amount of RAM. This is an empirical value observed for biggest mouse chromosomes.
    threads: 16 # Quick hack here to use a full node on Sacapus instead of just one core to avoid memory overload. Todo: replace that with a correct RAM resource definition.
    shell:"""
    #{input.bam2bedpe} {input.bam} > {output.bed}
    # Just replacing bam2bedpe with its 1-line content:
    {input.samtools} view -f 67 -F 1804 {input.bam} |  \
        awk '{{if ($9>0) print $3, $4, $4+$9, $5; if ($9<0) print $3, $8, $8-$9, $5}}' > {output.bed}

    {input.python} {input.puffin} {output.bed}
    """

rule samtools_view_get_chr_list:
    """
    Created: 2016-03-17 17h27
    Modified: 2016-12-08 18h46 - Adapted to new pattern and naming convention.
    
    Usage:
        "out/samtools/view/split_bam_by_chr/{filler}/SQ_list.txt"
    """
    input:
        samtools="opt/samtools-1.3.1/samtools",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        chr_list="out/samtools/view/get_chr_list/{filler}/SQ_list.txt"
    threads: 1
    shell:"""
    # Getting the header of the input
    # Retrieving the chromosomes in the header of the bam.
    # Contained in the 2nd column of the lines starting by @SQ. They start with "SN:" which is removed. 
    {input.samtools} view -H {input.bam} | grep @SQ | cut -f2 | sed 's/SN://g' > {output.chr_list}

    # Code below is now replaced by rule 'samtools_view_split_bam_by_chr'.
    #for CHR in `cat {output.chr_list}`
    #do
    #    {input.samtools} view -b {input.bam} $CHR > {output.chr_list}_$CHR.bam
    #    {input.samtools} index {output.chr_list}_$CHR.bam
    #done
    """

rule samtools_view_split_bam_by_chr:
    """
    Created: 2016-12-09 10h00 - chr list is now hardcoded in 'my_variables.snake' for easier handling in future steps, to avoid use of dynamics file outputs.
    """
    input:
        samtools="opt/samtools-1.3.1/samtools",
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai"
    output:
        bam="out/samtools/view/split_bam_by_chr/{filler}/{chr}.bam"
    threads: 1
    shell:"""
    {input.samtools} view -b {input.bam} {wildcards.chr} > {output.bam}
    """

