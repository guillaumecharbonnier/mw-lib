rule pysradb_download_gsm_se:
    """
    Created:
    Aim:
        It is more friendly to download SRA samples from GSM rather than (multiple) SRR ids to merge later.
        Pysradb is bugged right now for direct download of GSM (or I did not find the correct pipe of args)
        so I use pysradb to get the SRR from GSM, 
        then use fasterq-dump to get the fastq from SRR,
        and cat to merge all SRR files
    Test:
        out/pysradb_download_gsm_se/GSM2862905.fastq.gz
    """
    output:
        fastq="out/pysradb_download_gsm_se/{gsm}.fastq.gz"
    log:
        "out/pysradb_download_gsm_se/{gsm}.log"
    benchmark:
        "out/pysradb_download_gsm_se/{gsm}.benchmark.tsv"
    wildcard_constraints:
        gsm="GSM[0-9]+"
    conda:
        "../envs/pysradb_sra-tools.yaml"
    threads:
        1
    shell:
        """
        (
        SRR_IDS=`pysradb gsm-to-srr {wildcards.gsm} | cut -f2 -d " " | tail -n+2`
        DOWNLOAD_DIR="out/pysradb_download_gsm_se/{wildcards.gsm}"
        rm -rf $DOWNLOAD_DIR
        mkdir -p $DOWNLOAD_DIR
        cd $DOWNLOAD_DIR
        for ID in $SRR_IDS;
        do
            fasterq-dump $ID
        done
        gzip -c SRR*.fastq > ../{wildcards.gsm}.fastq.gz
        cd ../
        rm -rf {wildcards.gsm}
        ) &> {log}
        """


rule pysradb_download_gsm_pe:
    """
    Test:
        out/pysradb_download_gsm_pe/GSM5001892_1.fastq.gz
    """
    output:
        mates=expand("out/pysradb_download_gsm_pe/{{gsm}}_{mate}.fastq.gz", mate=["1","2"])
    log:
        "out/pysradb_download_gsm_pe/{gsm}.log"
    benchmark:
        "out/pysradb_download_gsm_pe/{gsm}.benchmark.tsv"
    wildcard_constraints:
        gsm="GSM[0-9]+"
    conda:
        "../envs/pysradb_sra-tools.yaml"
    threads:
        1
    shell:
        """
        (
        SRR_IDS=`pysradb gsm-to-srr {wildcards.gsm} | cut -f2 -d " " | tail -n+2`
        DOWNLOAD_DIR="out/pysradb_download_gsm_pe/{wildcards.gsm}"
        rm -rf $DOWNLOAD_DIR
        mkdir -p $DOWNLOAD_DIR
        cd $DOWNLOAD_DIR
        for ID in $SRR_IDS;
        do
            #fasterq-dump $ID
            fastq-dump $ID
        done
        gzip -c SRR*_1.fastq > ../{wildcards.gsm}_1.fastq.gz
        gzip -c SRR*_2.fastq > ../{wildcards.gsm}_2.fastq.gz
        cd ../
        rm -rf {wildcards.gsm}
        ) &> {log}
        """

#rule sra_tools_fastqdump_check_se_or_pe:
#    """
#    Created:
#        2016-09-30 9h40
#    Modified:
#        2017-02-14 11:20:18
#        2017-05-11 10:50:06 - Use miniconda sra-tools now.
#    Aim:
#        Print 1 if single-end and 2 if paired-end.
#        Sometimes better than looking at <LIBRARY_LAYOUT> in this kind of XML when experiment merge both designs.
#        https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRR1999050
#    Test:
#        out/sra-tools/fastq-dump_check_se_or_pe/SRR3938832.txt
#        out/sra-tools/fastq-dump_check_se_or_pe/SRR1202037.txt
#    """
#    output:
#        txt="out/sra-tools/fastq-dump_check_se_or_pe/{srr}.txt"
#    conda:
#        "../envs/sra-tools.yaml"
#    shell:
#        """
#        fastq-dump -I -X 1 -Z --split-spot {wildcards.srr} 2>/dev/null \
#            | awk '{{if(NR % 2 == 1) print substr($1,length($1),1)}}' \
#            | uniq \
#            | wc -l > {output.txt}
#        """
#
