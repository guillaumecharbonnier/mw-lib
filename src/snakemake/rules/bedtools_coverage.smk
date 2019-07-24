rule bedtools_coverage:
    """
    Created:
        2016-08-24 14h59
    Aim:
        Draft, maybe deprecated for my current need in comparison with deepTools multbigwigsummary + plotCorrelation
    """
    input:
        bam="out/bam/"
    output:
        ""
