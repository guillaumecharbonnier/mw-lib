rule bedtools_shuffle_various_tests:
    """
    Test:
        "out/bedtools/shuffle/various_tests_seed1/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions/bed_seed_incl.bed"
        
    """
    input:
        bed="out/{filler}.bed",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv",
        incl="out/awk/prepare_include_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        bed_seed_incl="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_seed_incl.bed",
        bed_seedpp_incl="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_seedpp_incl.bed",
        bed_seed_noincl="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_seed_noincl.bed",
        bed_seedpp_noincl="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_seedpp_noincl.bed",
        bed_noseed_incl_1="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_noseed_incl_1.bed",
        bed_noseed_incl_2="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_noseed_incl_2.bed",
        bed_noseed_noincl_1="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_noseed_noincl_1.bed",
        bed_noseed_noincl_2="out/bedtools/shuffle/various_tests_seed{seed}/{filler}/bed_noseed_noincl_2.bed",
    shell:"""
    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl} \
        -seed {wildcards.seed} > {output.bed_seed_incl}

    SEEDPP={wildcards.seed}+1
    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl} \
        -seed $SEEDPP > {output.bed_seedpp_incl}

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -seed {wildcards.seed} > {output.bed_seed_noincl}

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl} \
        -seed $SEEDPP > {output.bed_seedpp_noincl}

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl}  > {output.bed_noseed_incl_1}

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl}  > {output.bed_noseed_incl_2}

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} > {output.bed_noseed_noincl_1}

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} > {output.bed_noseed_noincl_2}

    """

rule bedtools_shuffle_noSeed:
    """
    Test:
        "out/bedtools/shuffle/noSeed1/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"
    """
    input:
        bed="out/{filler}.bed",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv",
        python="opt/miniconda/envs/py35/bin/python",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        bed="out/bedtools/shuffle/noSeed{noSeed}/{filler}.bed",
    wildcard_constraints:
        noSeed="[0-9]+"
    shell:"""
    # these two lines are very likely to be useless but as I have no idea how the seed is automatically taken, I just wait different time for each parallel jobs.
    NOSEED=`{input.python} -c "print(1/{wildcards.noSeed})"`
    sleep $NOSEED

    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} > {output.bed}
    """

rule bedtools_shuffle_noSeed_excl:
    """
    Test:
        "out/bedtools/shuffle/noSeed1_excl/bedtools/intersect/filter_danpos_ppr_on_mapability1/danpos/dtriple_v2/samtools/merge/samtools/sam_to_bam/bowtie2/pe_mm9/fastx-toolkit/fastx_trimmer/len30/merge_lanes/run113_run125_run126/PSK-SC-WT/danpos.smooth.positions.bed"
    """
    input:
        bed="out/{filler}.bed",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv",
        unmapable="out/bedtools/subtract/prepare_exclude_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bed",
        python="opt/miniconda/envs/py35/bin/python",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        bed="out/bedtools/shuffle/noSeed{noSeed}_excl/{filler}.bed",
    wildcard_constraints:
        noSeed="[0-9]+"
    shell:"""
    # these two lines are very likely to be useless but as I have no idea how the seed is automatically taken, I just wait different time for each parallel jobs.
    NOSEED=`{input.python} -c "print(1/{wildcards.noSeed})"`
    sleep $NOSEED

    {input.bedtools} shuffle \
        -i {input.bed} \
        -excl {input.unmapable} \
        -g {input.chromInfo} > {output.bed}
    """

rule bedtools_shuffle_seed:
    """
    Created: 2016-12-05 10h03
    First aim: shuffling danpos positions in order to get fraction of small structures sharing positions with nucleosomes.
    Note: limits for comparison with danpos positions: no exclusion of unmappable regions right now.
    + danpos maps in windows (50bp or more ?) whereas shuffle maps to 1bp window...

    Example:
        bed="out/bedtools/intersect/filter_danpos_{ppr}_on_mapability1/{filler}/{exp}_{sample}.bed"

    """
    input:
        bed="out/{filler}.bed",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv",
        incl="out/awk/prepare_include_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        bed="out/bedtools/shuffle/seed{seed}/{filler}.bed"
    threads: 1
    shell:"""
    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl} \
        -seed {wildcards.seed} > {output.bed}
    """

rule bedtools_shuffle_noOverlapping_seed:
    """
    Created: 2016-12-05 10h03
    First aim: shuffling danpos positions in order to get fraction of small structures sharing positions with nucleosomes.
    Note: limits for comparison with danpos positions: no exclusion of unmappable regions right now.
    + danpos maps in windows (50bp or more ?) whereas shuffle maps to 1bp window...

    Example:
        bed="out/bedtools/intersect/filter_danpos_{ppr}_on_mapability1/{filler}/{exp}_{sample}.bed"
    """
    input:
        bed="out/{filler}.bed",
        chromInfo="out/mysql/ucsc/mm9/chromInfo_main.tsv",
        incl="out/awk/prepare_include_for_shuffle/ucsc/bigWigToBedGraph/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bedGraph",
        bedtools="opt/bedtools2/bin/bedtools"
    output:
        bed="out/bedtools/shuffle/noOverlapping_seed{seed}/{filler}.bed"
    threads: 1
    shell:"""
    {input.bedtools} shuffle \
        -i {input.bed} \
        -g {input.chromInfo} \
        -incl {input.incl} \
        -seed {wildcards.seed} \
        -noOverlapping > {output.bed}
    """

