ruleorder: rsync_ucsc_cse > rsync_extra

rule rsync_ucsc_cse:
    """
    Created:
        2017-01
    Aim:
        UCSC recommends the use of rsync to retrieve data from their FTP server.
    Note:
        I think 'cse' is one mirror. But there are others like 'soe'.
    Test:
        out/rsync/ucsc/gbdb/mm9/bbi/crgMapabilityAlign75mer.bw
        out/rsync/ucsc/gbdb/mm9/liftOver/mm9ToMm10.over.chain.gz
        out/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt.gz
        out/rsync/ucsc/goldenPath/mm10/database/cpgIslandExt.txt.gz
        out/rsync/ucsc/goldenPath/mm10/database/rmsk.txt.gz
        out/rsync/ucsc/goldenPath/mm10/bigZips/mm10.2bit
        out/rsync/ucsc/goldenPath/mm10/database/cytoBand.txt.gz
    """
    output:
        "out/rsync/ucsc/{ucsc_ftp_file}"
    params:
        mirror="cse" # You can try 'soe'
    shell:
        "rsync -a -P rsync://hgdownload.cse.ucsc.edu/{wildcards.ucsc_ftp_file} {output}"

rule rsync_extra:
    """
    Created:
        2019-01-30 11:17:36
    Test:
        out/rsync/ucsc/goldenPath/hg19/database/tfbsConsSites.txt.gz
        out/rsync/_-aP/hgdownload.cse.ucsc.edu/goldenPath/hg19/database/tfbsConsSites.txt.gz
        out/rsync/_-aP/hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit
    """
    output:
        "out/{tool}{extra}/{url}"
    wildcard_constraints:
        tool="rsync/"
    params:
        extra = params_extra
    shell:
        "rsync {params.extra} rsync://{wildcards.url} {output}"

