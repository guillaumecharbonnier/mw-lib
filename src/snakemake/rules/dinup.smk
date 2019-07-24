rule dinup:
    """
    Created:
        2016-03-17
    Modified:
        2016-03-21 - From test to rule usable in parallel.
    TODO:
        Test this rule
    """
    input:
        bed1="out/bedtools/bamtobed/{index}/{exp}/{sample1}.bed",
        bed2="out/bedtools/bamtobed/{index}/{exp}/{sample2}.bed"
    output:
        done="dinup_{index}_{exp}_{sample1}_{sample2}.done"
    threads:
        1
    shell:
        """
        WDPATH=`pwd`
        NAME="{wildcards.index}_{wildcards.exp}_{wildcards.sample1}_{wildcards.sample2}"
        TMPPATH="tmp/dinup/$NAME"
        
        mkdir --parents $TMPPATH
        cd $TMPPATH
        python2.7 $WDPATH/opt/dinup_1.3/bin/dinup -t $WDPATH/{input.bed1} -c $WDPATH/{input.bed2} -n $NAME
    
        mv $TMPPATH/DiNuP_results/$NAME_dinup.bed out/dinup/{wildcards.index}/{wildcards.exp}/{wildcards.sample1}_{wildcards.sample2}/dinup.bed
        mv $TMPPATH/DiNuP_results/$NAME_dinup.xls out/dinup/{wildcards.index}/{wildcards.exp}/{wildcards.sample1}_{wildcards.sample2}/dinup.xls
        
        touch $WDPATH/{output.done}
        """

