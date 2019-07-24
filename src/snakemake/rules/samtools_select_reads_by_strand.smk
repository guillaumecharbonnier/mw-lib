rule select_reads_by_strand:
    """
    Created: 2016-11-14 15h31 - From Denis' workflow.
    """
    input:
        bam="out/{id}.bam",
    output:
        min="out/samtools/select_reads_by_strand/{id}_min.bam",
        plus="out/samtools/select_reads_by_strand/{id}_plus.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell: 
        """
        samtools view -f99 -hb {input.bam} > {output.min}_99.bam
        samtools view -f147 -hb {input.bam} > {output.min}_147.bam
        samtools view -f83 -hb {input.bam} > {output.plus}_83.bam
        samtools view -f163 -hb {input.bam} > {output.plus}_163.bam

        samtools merge -f -h {output.min}_99.bam {output.min} {output.min}_99.bam {output.min}_147.bam
        samtools merge -f -h {output.plus}_83.bam {output.plus} {output.plus}_83.bam {output.plus}_163.bam

        rm -f {output.min}_99.bam {output.min}_147.bam {output.plus}_83.bam {output.plus}_163.bam
        """
