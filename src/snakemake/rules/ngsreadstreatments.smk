rule ngsreadstreatment_se:
    """
    Test:
        out/ngsreadstreatment/se/sra-tools/fastq-dump_se/SRR3126243.fastq
    """
    input:
        jar="out/wget/https/sourceforge.net/projects/ngsreadstreatment/files/NgsReadsTreatment_v1.1.jar",
        fastq="out/{filler}.fastq"
    output:
        fastq="out/ngsreadstreatment/se/{filler}.fastq",
        report="out/ngsreadstreatment/se/{filler}.txt"
    log:
        "out/ngsreadstreatment/se/{filler}.log"
    shell:
        """
        (mkdir -p {output.fastq}_tmp
        ln -f {input.fastq} {output.fastq}_tmp/sample.fastq
        java -jar {input.jar} {output.fastq}_tmp/sample.fastq
        mv {output.fastq}_tmp/sample_1_trated.fastq {output.fastq}
        mv {output.fastq}_tmp/Report.log {output.report} ) &> {log}
        """

