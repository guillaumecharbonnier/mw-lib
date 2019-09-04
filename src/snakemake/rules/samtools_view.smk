rule samtools_view_sam_to_bam_extra:
    """
    Created:
        2019-01-31 22:44:26
    Aim:
        Convert sam to bam. Additionnal filters can be applied using params_extra.
    Test:
        out/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
        out/samtools/view_sam_to_bam_-F_4_-q_5/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
    """
    input:
        sam = "out/{filler}.sam"
    output:
        bam = "out/{tool}{extra}/{filler}.bam"
    log:
        "out/{tool}{extra}/{filler}.log"
    wildcard_constraints:
        tool="samtools/view_sam_to_bam"
    params:
        extra = params_extra
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -bSh {params.extra} {input.sam} -o {output.bam} &> {log}"

rule samtools_view_bam_to_bam_extra:
    """
    Created:
        2019-02-10 20:04:32
    Aim:
        Allow to apply filter on bam files.
    Test:
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
    output:
        bam = "out/{tool}{extra}/{filler}.bam"
    log:
        "out/{tool}{extra}/{filler}.log"
    wildcard_constraints:
        tool="samtools/view_bam_to_bam"
    params:
        extra = params_extra
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -bh {params.extra} {input.bam} -o {output.bam} &> {log}"

rule samtools_view_bam_to_sam_extra:
    """
    Created:
        2019-02-10 20:04:32
    Aim:
        Allow to apply filter on bam files.
    Test:
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
    output:
        sam = "out/{tool}{extra}/{filler}.sam"
    log:
        "out/{tool}{extra}/{filler}.log"
    wildcard_constraints:
        tool="samtools/view_bam_to_sam"
    params:
        extra = params_extra
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -h {params.extra} {input.bam} -o {output.sam} &> {log}"


rule samtools_view_bam_to_txt_extra:
    """
    Created:
        2019-02-10 20:04:32
    Aim:
        For some combinations of arguments, samtools view outputs text file.
    Test:
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
    output:
        txt = "out/{tool}{extra}/{filler}.txt"
    log:
        "out/{tool}{extra}/{filler}.log"
    wildcard_constraints:
        tool="samtools/view_bam_to_txt"
    params:
        extra = params_extra
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view {params.extra} {input.bam} -o {output.txt} &> {log}"

rule samtools_view_filter_bam_with_bed:
    """
    Created:
        2018-04-04 11:23:07
    Aim:
        Extracting only autosomes and sexual chromsomes.
    Test:
        out/samtools/view_filter_bam_with_bed_bed-mm10-main-chr/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243.bam
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
        bed = lambda wildcards: config['ids'][wildcards.bed_id],
    output:
        bam="out/{tool}{extra}_{bed_id}/{filler}.bam"
    wildcard_constraints:
        tool="samtools/view_filter_bam_with_bed"
    params:
        extra = params_extra
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools view -bh {params.extra} -L {input.bed} {input.bam} > {output.bam}"




#Legacy past this point.
#Keeping them here for now because I have not tried to implement them using extra rules above.


# Generalization of samtools view is not needed right now
#rule samtools_view_extra_sam_to_bam:
#    """
#    Created:
#        2018-12-06 11:12:17
#    Aim:
#        Convert sam to bam and optionally apply filtering
#    Common uses:
#        Converting sam to bam while preserving header:
#            out/samtools/view_-bSh/{filler}.bam
#        Converting sam to bam while preserving header and filter on quality:
#            out/samtools/view_-bSh_-q_20/{filler}.bam
#    Doc:
#        samtools view [options] in.bam|in.sam|in.cram [region...]
#
#        With no options or regions specified, prints all alignments in the specified input alignment file (in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).
#
#        You may specify one or more space-separated region specifications after the input filename to restrict output to only those alignments which overlap the specified region(s). Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format).
#
#        The -b, -C, -1, -u, -h, -H, and -c options change the output format from the default of headerless SAM, and the -o and -U options set the output file name(s).
#
#        The -t and -T options provide additional reference data. One of these two options is required when SAM input does not contain @SQ headers, and the -T option is required whenever writing CRAM output.
#
#        The -L, -r, -R, -q, -l, -m, -f, and -F options filter the alignments that will be included in the output to only those alignments that match certain criteria.
#
#        The -x, -B, and -s options modify the data which is contained in each alignment.
#
#        Finally, the -@ option can be used to allocate additional threads to be used for compression, and the -? option requests a long help message.
#    Tests:
#        out/samtools/view_sam_to_bam_-bSh/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_pe_raw/inp/fastq/run178/S001797_BAC_mT_DHS_PGK_p5424_rep2-49056/BAC-mT-DHS-PGK-p5424-rep2_S8.bam
#    """
#    input:
#        sam = "out/{filler}.{extin}"
#    output:
#        #bam="out/{tool}{extra}/{filler}.bam"
#        bam="out/{tool}_{extin}_to_{extout}_{extra}{L_bedid}/{filler}.{extout}"
#    wildcard_constraints:
#        tool="samtools/view",
#        extin='sam|bam',
#        extout='sam|bam',
#    params:
#        extra = params_extra
#    conda:
#        "../envs/samtools.yaml"
#    shell:
#        "samtools view {params.extra} {input.bed} {input.sam_or_bam} > {output.sam_or_bam}"


# Attempt to generalize samtools view for all usecase but it happens to be messier than having just a few different rules for each usecase.
#def input_bai_samtools_view(wildcards):
#    if wildcards['extin'] == 'bam':
#        path = 'out/' + wildcards['filler'] + '.bam.bai'
#        return(path)
#    else:
#        # Returning input because snakemake is not happy if I return nothing
#        path ='out/' + wildcards['filler'] + '.' + wildcards['extin']
#    return(path)
#
#def input_bed_samtools_view(wildcards):
#    if wildcards['L_bedid'] == '':
#        # Returning input because snakemake is not happy if I return nothing
#        path = 'out/' + wildcards['filler'] + '.' + wildcards['extin']
#    else:
#        print('write code in input_bed_samtools_view')
#        # strip _-L_ then call global function input_bed
#        # to return path to bed file.
#    return(path)
#
# TODO: Write params_L_bedpath function here
#
#rule samtools_view_extra_sam_to_bam:
#    """
#    Created:
#        2018-12-06 11:12:17
#    Aim:
#        Convert sam to bam and optionally apply filtering
#    Common uses:
#        Converting sam to bam while preserving header:
#            out/samtools/view_-bSh/{filler}.bam
#        Converting sam to bam while preserving header and filter on quality:
#            out/samtools/view_-bSh_-q_20/{filler}.bam
#    Tests:
#        out/samtools/view_sam_to_bam_-bSh/bowtie2/pe_GRCm38/sickle/pe_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_pe_raw/inp/fastq/run178/S001797_BAC_mT_DHS_PGK_p5424_rep2-49056/BAC-mT-DHS-PGK-p5424-rep2_S8.bam
#    """
#    input:
#        #sam = "out/{filler}.sam"#,
#        sam_or_bam = "out/{filler}.{extin}",
#        bai = input_bai_samtools_view, #Return {input.sam_or_bam} if {extin} is sam.
#        bed = input_bed_samtools_view  #Return {input.sam_or_bam} if {L_bedid} is empty.
#    output:
#        #bam="out/{tool}{extra}/{filler}.bam"
#        sam_or_bam="out/{tool}_{extin}_to_{extout}_{extra}{L_bedid}/{filler}.{extout}"
#    wildcard_constraints:
#        tool="samtools/view",
#        extin='sam|bam',
#        extout='sam|bam',
#        extra='-*[^(\/|_\-L_)]*', # extra could be empty or starting with '-' and not matching '/' (begining of {filler} or '_-L_' (begining of {L_bedid}
#        L_bedid='|_-L_[a-zA-Z0-9-]+' # L_bedid could be empty or starting with '_-L_' then 'bedid'.
#    params:
#        extra = params_extra,
#        L_bedpath = params_L_bedpath
#    conda:
#        "../envs/samtools.yaml"
#    shell:
#        "samtools view {params.extra} {params.L_bedpath} {input.sam_or_bam} > {output.sam_or_bam}"

rule samtools_view_bSh:
    """
    Created:
        2017-05-06 22:50:39
    Aim:
        Convert sam to bam without sorting.
    Recommended follow-up:
        samtools_sort
    Test:
        out/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR1202037.bam
    """
    input:
        sam="out/{filler}.sam"
    output:
        bam="out/samtools/view_bSh/{filler}.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        "samtools view -bSh {input.sam} > {output.bam}"

rule samtools_view_bSh_q_legacy:
    """
    Created:
        2017-10-27 17:46:06
    Aim:
        Convert sam to bam without sorting and filter on quality.
    Recommended follow-up:
        samtools_sort
    """
    input:
        sam="out/{filler}.sam"
    output:
        bam="out/samtools/view_bSh_q-{q}/{filler}.bam"
    wildcard_constraints:
        q="[0-9]+"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    shell:
        """
        samtools view -bSh -q {wildcards.q} {input.sam} > {output.bam}
        """

rule samtools_view_b_F_q:
    """
    Created:
        2018-03-16 18:00:34
    Aim:
        Filtering according to Blueprint protocol.
        http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38
    Test:
        out/samtools/view_b_F-4_q-5/picard/MarkDuplicates_RemoveDuplicates-false_AssumeSorted-true/picard/SortSam_sortOrder-coordinate/samtools/view_bSh/bwa/samse_q-5_fa-GRCh38-Blueprint/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.bam
    """
    input:
        bam="out/{filler}.bam"
    output:
        bam="out/samtools/view_b_F-{F}_q-{q}/{filler}.bam"
    wildcard_constraints:
        q="[0-9]+",
        F="[0-9]+"
    threads:
        1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools\
            view\
            -b\
            -F {wildcards.F}\
            -q {wildcards.q}\
            {input.bam} > {output.bam}
        """

rule samtools_view_b_F:
    """
    Created:
        2018-03-16 18:00:34
    Aim:
        Filtering according to Blueprint protocol.
        http://dcc.blueprint-epigenome.eu/#/md/chip_seq_grch38
    Test:
        out/samtools/view_b_F-1024/samtools/view_b_F-4_q-5/picard/MarkDuplicates_RemoveDuplicates-false_AssumeSorted-true/picard/SortSam_sortOrder-coordinate/samtools/view_bSh/bwa/samse_q-5_fa-GRCh38-Blueprint/gunzip/merge_lanes_nextseq500_se_raw/inp/fastq/run170/S001580_TH125_EC_K27ac-41046/TH125-EC-K27ac_S6.bam
    """
    input:
        bam="out/{filler}.bam"
    output:
        bam="out/samtools/view_b_F-{F}/{filler}.bam"
    wildcard_constraints:
        F="[0-9]+"
    threads:
        1
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools\
            view\
            -b\
            -F {wildcards.F}\
            {input.bam} > {output.bam}
        """



rule samtools_view_bh_chr:
    """
    Created:
        2018-04-04 11:23:07
    Aim:
        Extracting only one chromsome
    Test:
        out/gemBS/bscall_blueprint_like/samtools/view_bh_chr19/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028/done

    """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        samtools="opt/miniconda/envs/samtools/bin/samtools"
    output:
        bam="out/samtools/view_bh_{chr}/{filler}.bam"
    wildcard_constraints:
        chr="[a-zA-Z0-9]+"
    threads:
        1
    shell:
        """
        {input.samtools}\
            view\
            -bh\
            {input.bam}\
            {wildcards.chr} > {output.bam}
        """


rule samtools_view_bh_L:
    """
    Created:
        2018-04-04 11:23:07
    Aim:
        Extracting only autosomes and sexual chromsomes.
        Done because I want to get rid of chrEBV in hg38 alignemnts.
    Note:
        samtools view -h aln.bam | awk '{if($3 != "chrX" && $3 != "chrY"){print $0}}' | samtools view -Sb - > aln.filter.bam

        A likely faster method might be to just make a BED file containing those chromosomes/contigs and then just:
        samtools view -b -L chromosomes.bed foo.bam > subset.bam

    Test:
        out/samtools/view_bh_L-hg38-main-chr/samtools/sort/samtools/view_bSh_q-30/bowtie2/pe_GRCh38_--very-sensitive-local/sickle/pe_-t_sanger_-q_20/gunzip/merge_lanes_nextseq500_pe_raw/inp/fastq/run223/S002289_H4K5ac_REH_gev-91095/H4K5ac-REH-gev_S5.bam
    """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
        bed = lambda wildcards: config['ids'][wildcards.bed_id],
        samtools="opt/miniconda/envs/samtools/bin/samtools"
    output:
        bam="out/samtools/view_bh_L-{bed_id}/{filler}.bam"
    threads:
        1
    shell:
        """
        {input.samtools} view -bh -L {input.bed} {input.bam} > {output.bam}
        #{input.samtools} view -h {input.bam} |\
        #    awk '{{ if($3 ~ "^chr[0-9XY]+$"){{print $0}} }}' |\
        #        {input.samtools} view -bSh - > {output.bam}
        """


rule samtools_filterQ30:
        """
        Filter aligned reads on mapping quality criterium.

	Usage:
/bedtools/intersect/bam_in_feature/mm9/run113/
	expand("data/samtools_filterQ30/{id}_Q30.bam", id="bedtools/intersect/bam_in_feature/mm9/run113/1-MNS-TGC-WT/simpleRepeats")

data/samtools_filterQ30/bedtools/intersect/bam_in_feature/mm9/run113/1-MNS-TGC-WT/simpleRepeats_Q30.bam
        """
        input:
            "out/{id}.bam"
        output:
            "out/samtools/filterQ30/{id}_Q30.bam"
        shell:
            """
            samtools view -bh -q 30 {input} > {output}
            """


rule samtools_view_check_se_or_pe:
    """
    Created:
        2018-01-15 17:37:43
    Aim:
        Check if a bam file is single-end (output 0) or paired-end (output 1)
    Note:
        samtools view -c -f 1
    """
    input:
        samtools="opt/miniconda/envs/samtools/bin/samtools",
        bam="out/{filler}.bam"
    output:
        txt="out/samtools/view_check_se_or_pe/{filler}.txt",
    shell:
        """
        {input.samtools} view -c -f 1 {input.bam} > {output.txt
        """

rule samtools_view_paired_end_to_single_end:
    """
    Created:
        2016-07-18 18h13
    Aim:
        Select only the reads that are first in a pair.
        Bam produced contain reads and consistent size but nothing is displayed in IGV. I am not sure macs1.4 will take them as single-end for peak-calling so I will compare with macs2 paired-end peak-calling.
    """
    input:
        samtools="opt/samtools-1.3.1/samtools",
        bam="out/bam/{index}/{run}/{sample}.bam"
    output:
        bam="out/samtools/pe_to_se/{index}/{run}/{sample}.bam",
        ln="out/bam/{index}/{run}_pe_to_se/{sample}.bam"
    shell:
        """
        {input.samtools} view -bSh -f 0x0040 {input.bam} | {input.samtools} sort -@ {threads} -m 10G - > {output.bam}
        {input.samtools} index {output.bam}
        ln {output.bam} {output.ln}
        ln {output.bam}.bai {output.ln}.bai
        """
"""
Modified:
    2017-05-10 18:18:14 - These rules are here for legacy. I recommend using 'samtools/sort/samtools/view_bSh' instead of 'samtools/sam_to_bam'. It is faster and more respectful of the coding conventions.
"""

rule samtools_sam_to_bam:
    """
    Created:
        2016-09-28 11:40 - Added when I had to analyze huge MNase data from GEO.
    Modified:
        2016-09-30 11h53 - Added prefix with T argument to avoid flooding my working directory. Not tested.
        2017-02-20 15:30:15 - Modified input/output patterns.
        2017-05-04 16:31:32 - Reduced sam size limit to use pipe mode from 100000000 to 10000000
    Note:
        Pipe are not suitable for big files and produce broken pipe errors.
        I set the threshold at 100go but this was empiric.
        Actually I am not even sure that piping is faster for standard files.
    """
    input:
        sam="out/{filler}.sam",
        samtools="opt/samtools-1.3.1/samtools"
    output:
        bam="out/samtools/sam_to_bam/{filler}.bam"
    threads: 1
    shell:
        """
        SAM_SIZE=`du {input.sam} | cut -f1`

        if (($SAM_SIZE < 10000000));
        then
            echo "Pipe mode."
            {input.samtools} view -bSh {input.sam} | \
            {input.samtools} sort -@ {threads} -T {output.bam} -m 10G - > {output.bam}
        else
            echo "Huge file detected. Temporary file mode."
            {input.samtools} view -bSh {input.sam} > {output.bam}_unsorted.bam
            {input.samtools} sort -@ {threads} -T {output.bam} -m 10G {output.bam}_unsorted.bam > {output.bam}
            rm -f {output.bam}_unsorted.bam
        fi
        """

rule samtools_sam_to_bam_legacy:
    """
    Created: 2016-02-04 11:40
    Modified: 2016-09-30 11h53 - Added prefix with T argument to avoid flooding my working directory. Not tested.
    Modified: 2016-11-24 10h26 - Changed to 'legacy' because the new pattern concept makes this rule obsolete.
    """
    input:
        sam="out/sam/{index}/{run}/{sample}.sam",
        samtools="opt/samtools-1.3.1/samtools"
    output:
        bam="out/bam/{index}/{run}/{sample}.bam"
    threads: 1
    shell:
        """
        {input.samtools} view -bSh {input.sam} | \
            {input.samtools} sort -@ {threads} -T {output.bam} -m 10G - > {output.bam}
        """

rule samtools_view_b_h_s:
    """
    Created:
        2018-05-28 11:11:02
    Aim:
        Subsample bam
    Note:
        -s {wildcards.seed}.${{SCALEFACTOR}}
    Test:
        out/samtools/view_b_h_s-0.5/ln/alias/experiments/hg38_H3K27ac_thymus/CD34.bam
        out/samtools/view_b_h_s-1.5/ln/alias/experiments/hg38_H3K27ac_thymus/CD34.bam

    """
    input:
        samtools="opt/miniconda/envs/samtools/bin/samtools",
        bam="out/{filler}.bam"
    output:
        bam="out/samtools/view_b_h_s-{s}/{filler}.bam"
    wildcard_constraints:
        s="[0-9].[0-9]+"
    threads:
        16
    shell:
        """
        {input.samtools} view -bh -s {wildcards.s} --threads {threads} {input.bam} > {output.bam}
        """


rule samtools_view_subsampling_bam_get_sample_size:
    """
    Created: 2016-06-03 9h34
    The aim is to assess the number of mapped read in multiple condition then do subsampling in the conditions to reach a given number of read for all samples and multiple replicates for each sample.

    The rule is splitted into two parts "get_scale_factor_to_100m" and "produce_subsamples" to allow parallel subsampling with different seeds "-s" without redoing each time the scale factor computation which is quite time consuming per se.
    """
    input:
        bam="out/{id}.bam",
        samtools="opt/samtools-1.3.1/samtools"
    output:
        txt="out/samtools/view/subsampling_bam/n_paired_fragments/{id}.txt"
    shell:"""
    {input.samtools} view -c -f 0x2 -F 0x4 -q 10 {input.bam} > {output.txt}
    # Here it may be better to use '-f 0x2' instead of '-F 0x40' if data is paired end (see discussion in rule below).
    # I have to check what is selected using '-f 0x2' if data is single end.
    """

rule samtools_view_subsampling_bam_produce_subsamples:
    """
    Created: 2016-06-06 15h54
    Warning: output is single-end even if input is paired-end. I am not sure if the output file is biased depending on the way the two paired-end mates could be picked. Independent picking? Ignored for the moment.
    """
    input:
        bam="out/{id}.bam",
        txt="out/samtools/view/subsampling_bam/n_paired_fragments/{id}.txt",
        samtools="opt/samtools-1.3.1/samtools"
    output:
        bam="out/samtools/view/subsampling_bam/seed{seed}_{subsamplesize}pf/{id}.bam",
        bai="out/samtools/view/subsampling_bam/seed{seed}_{subsamplesize}pf/{id}.bam.bai"
    threads: 2
    wildcard_constraints:
        seed="[0-9]+",
        subsamplesize="[0-9]+"
    shell:"""
    # The scale factor is actually trimmed with "cut" to get only the decimal part as required by samtools.
    SCALEFACTOR=`awk '{{print {wildcards.subsamplesize}/$1}}' {input.txt} | cut -f2 -d "."`

    # The flag filtering is supposed to remove unmapped reads.
    #{input.samtools} view -bh -F 0x40 --threads {threads} -s {wildcards.seed}.${{SCALEFACTOR}} {input.bam} > {output.bam}

    # Suspecting the -F 0x40 to mess up with paired end selection.
    #{input.samtools} view -bh --threads {threads} -s {wildcards.seed}.${{SCALEFACTOR}} {input.bam} > {output.bam}

    # As the previous testing revealed the issue with mixing -s and -F, I try to do the two steps sequentially instead of in the same command:
    # This has to be test, e.g with seed3 patterns.
    #{input.samtools} view -bh -F 0x40 --threads {threads} {input.bam} | \
    #    {input.samtools} view -bh -s {wildcards.seed}.${{SCALEFACTOR}} - > {output.bam}

    #Apparently, the command above still produce single-end reads. Maybe I should try other flags for filtering.
    # Testing with seed4 Input WT
    #{input.samtools} view -bh -F 4 --threads {threads} {input.bam} | \
    #    {input.samtools} view -bh -s {wildcards.seed}.${{SCALEFACTOR}} - > {output.bam}

    # The previous test is quite succesful as I am able to keep paired-end reads. But I also find some read with unmapped mate. It may be interesting to only subsample on correctly paired reads...
    # The -f 0x2 part will get only "properly paired" alignments:
    # Testing with seed4 H5K5acWT
    #{input.samtools} view -bh -f 0x2 --threads {threads} {input.bam} | \
    #    {input.samtools} view -bh -s {wildcards.seed}.${{SCALEFACTOR}} - > {output.bam}

    # The last test is not successful, because i can still see some reads with unmapped pairs (or interchromosomal pairs, but this is expected).
    # If the last test is successful, it remains to test if the samtools view can be one line: does it select proper paired reads then apply the subsampling (wanted behavior) or does these steps in the opposite order (unwanted).
    # I try to add both '-f 0x2' and '-F 0x4' according to this post:
    # https://www.biostars.org/p/115437/ (Istvan Albert's post)
    #{input.samtools} view -bh -f 0x2 -F 0x4 -q 10 -s {wildcards.seed}.${{SCALEFACTOR}} --threads {threads} {input.bam} > {output.bam}

    # The previous line (with adjusted get_sample_size rule) produces approx 2200000 fragment when I require 5000000. Now I check if I get closer result if I split this in two command:
    # Tested with seed 7:
    #{input.samtools} view -bh -f 0x2 -F 0x4 -q 10 --threads {threads} {input.bam} | \
    #    {input.samtools} view -bh -s {wildcards.seed}.${{SCALEFACTOR}} - > {output.bam}

    # The previous line produces also approx 2200000 fragment when I require 5000000. Thus the two-line command does not solve the issue. This is likely to be an error in the scale factor calculation. Maybe it counts reads instead of fragments? Anyway it does not affect downward analysis as the normalization by subsampling is correctly done.
    {input.samtools} view -bh -f 0x2 -F 0x4 -q 10 -s {wildcards.seed}.${{SCALEFACTOR}} --threads {threads} {input.bam} > {output.bam}

    {input.samtools} index {output.bam} {output.bai}
    """

# Below are deprecated rules

rule samtools_view_subsampling_bam_get_sample_size_old:
    """
    Created: 2016-06-03 9h34
    The aim is to assess the number of mapped read in multiple condition then do subsampling in the conditions to reach a given number of read for all samples and multiple replicates for each sample.

    The rule is splitted into two parts "get_scale_factor_to_100m" and "produce_subsamples" to allow parallel subsampling with different seeds "-s" without redoing each time the scale factor computation which is quite time consuming per se.
    """
    input:
        bam="out/bam/{index}/{exp}/{sample}.bam",
        samtools="opt/samtools-1.3.1/samtools"
    output:
        txt="out/samtools/view/subsampling_bam/{index}/{exp}/{sample}/n_paired_fragments.txt"
    shell:"""
    {input.samtools} view -c -F 0x40 {input.bam} > {output.txt}
    """

rule samtools_view_subsampling_bam_produce_subsamples_old:
    """
    Created: 2016-06-06 15h54

    Remark about threads: When I give 4 threads, it uses only 200% CPU so I give him only 2 threads.
    """
    input:
        bam="out/bam/{index}/{exp}/{sample}.bam",
        txt="out/samtools/view/subsampling_bam/{index}/{exp}/{sample}/n_paired_fragments.txt",
        samtools="opt/samtools-1.3.1/samtools"
    output:
        bam="out/samtools/view/subsampling_bam/{index}/{exp}/{sample}/seed{seed,[0-9]+}_{subsamplesize, [0-9]+}pf.bam",
        bai="out/samtools/view/subsampling_bam/{index}/{exp}/{sample}/seed{seed}_{subsamplesize}pf.bam.bai",
        ln_bam="out/bam/{index}/subsampling_{exp}/{sample}_s{seed}_{subsamplesize}pf.bam",
        ln_bai="out/bam/{index}/subsampling_{exp}/{sample}_s{seed}_{subsamplesize}pf.bam.bai"
    threads: 2
    shell:"""
    # The scale factor is actually trimmed with "cut" to get only the decimal part as required by samtools.
    SCALEFACTOR=`awk '{{print {wildcards.subsamplesize}/$1}}' {input.txt} | cut -f2 -d "."`
    {input.samtools} view -b -F 0x40 --threads {threads} -s {wildcards.seed}.${{SCALEFACTOR}} {input.bam} > {output.bam}
    {input.samtools} index {output.bam} {output.bai}

    ln {output.bam} {output.ln_bam}
    ln {output.bai} {output.ln_bai}
    """
