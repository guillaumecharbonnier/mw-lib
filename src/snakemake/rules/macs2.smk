rule macs2_callpeak_extra:
    """
    Created:
        2018-10-22 12:33:01
    Aim:
        Call peaks with control/input
    Doc:
        https://github.com/taoliu/MACS
    Test:
        out/macs2/callpeak_--broad/samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/bowtie2/se_mm10/sickle/se_-t_sanger_-q_30/sra-tools/fastq-dump_se/SRR3126243_over_SRR3126242_peaks.bed
    """
    input:
        bam_chip  = "out/{filler}/{chip}.bam",
        bai_chip  = "out/{filler}/{chip}.bam.bai",
        bam_input = "out/{filler}/{input}.bam",
        bai_input = "out/{filler}/{input}.bam.bai"
    output:
        bed    = "out/{tool}{extra}/{filler}/{chip}_over_{input}_peaks.bed",
        xls    = "out/{tool}{extra}/{filler}/{chip}_over_{input}_peaks.xls"
    log:
                 "out/{tool}{extra}/{filler}/{chip}_over_{input}_peaks.log"
    benchmark:
                 "out/{tool}{extra}/{filler}/{chip}_over_{input}_peaks.benchmark.tsv"
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra  = params_extra
    wildcard_constraints:
        tool = "macs2/callpeak"
    conda:
        "../envs/macs2.yaml"
    wildcard_constraints:
        chip  = "[a-zA-Z0-9-_]+",
        input = "[a-zA-Z0-9-_]+"
    shell:
        """
        (macs2 callpeak {params.extra}\
            --treatment {input.bam_chip}\
            --control {input.bam_input}\
            --name {wildcards.chip}_over_{wildcards.input}\
            --outdir {params.outdir}
        # Renaming narrowPeak or broadPeak to have only one output bed name for all variations of settings
        TO_RENAME=`find {params.outdir} -name '{wildcards.chip}_over_{wildcards.input}_peaks.narrowPeak' -o -name '{wildcards.chip}_over_{wildcards.input}_peaks.broadPeak'`
        echo $TO_RENAME
        # We want the output bed to have standard columns
        # The score column is limited to 1000
        awk 'BEGIN {{OFS = FS = "\\t"}} {{print $1,$2,$3,$4,$5,($6>1000)? 1000 : $6}}' $TO_RENAME > {output.bed}
        rm -f $TO_RENAME) &> {log}
        """

rule macs2_noctrl_callpeak_extra:
    """
    Created:
        2018-10-22 12:33:01
    Aim:
        Call peaks without control/input.
    Test:
        out/macs2/noctrl_callpeak_--format_BAM_--gsize_mm/samtools/sort/samtools/view_bSh/bowtie2/se_mm10/sickle/se_-t_sanger_-q_20/sra-tools/fastq-dump/SRR1202037_peaks.bed
    """
    input:
        bam_chip  = "out/{filler}/{chip}.bam",
        bai_chip  = "out/{filler}/{chip}.bam.bai",
    output:
        bed    = "out/{tool}{extra}/{filler}/{chip}_peaks.bed",
        xls    = "out/{tool}{extra}/{filler}/{chip}_peaks.xls"
    log:
                 "out/{tool}{extra}/{filler}/{chip}_peaks.log"
    benchmark:
                 "out/{tool}{extra}/{filler}/{chip}_peaks.benchmark.tsv"
    params:
        outdir = "out/{tool}{extra}/{filler}",
        extra  = params_extra
    wildcard_constraints:
        tool = "macs2/noctrl_callpeak"
    conda:
        "../envs/macs2.yaml"
    wildcard_constraints:
        chip="[a-zA-Z0-9-_]+",
        input="[a-zA-Z0-9-_]+"
    shell:
        """
        (macs2 callpeak {params.extra}\
            --treatment {input.bam_chip}\
            --name {wildcards.chip}\
            --outdir {params.outdir}
        # Renaming narrowPeak or broadPeak to have only one output bed name for all variations of settings
        TO_RENAME=`find {params.outdir} -name '{wildcards.chip}_peaks.narrowPeak' -o -name '{wildcards.chip}_peaks.broadPeak'`
        # We want the output bed to have standard columns
        # The score column is limited to 1000
        awk 'BEGIN {{OFS = FS = "\\t"}} {{print $1,$2,$3,$4,$5,($6>1000)? 1000 : $6}}' $TO_RENAME > {output.bed}
        rm -f $TO_RENAME) &> {log}
        """

# Legacy rules past this line.

rule macs2_callpeak_broad_format_gsize_keepdup_no_control_legacy:
    """
    Created:
        2018-02-10 00:14:58
    Aim:
        Straightforward peak calling, but not control by input...
        For ATAC:
        macs2 callpeak --treatment my.bam -g hs --bdg -q 0.05 --keep-dup all -f BAMPE --broad
    Test:
        out/macs2/callpeak_broad_format-BAM_gsize-hs_keepdup-all_no_control/ln/alias/experiments/hg38_ATAC_thymus/CD34_peaks.broadPeak
    """
    input:
        macs = "opt/miniconda/envs/macs2/bin/macs2",
        bam =  "out/{filler}/{sample}.bam",
        bai =  "out/{filler}/{sample}.bam.bai"
    output:
        bed =  "out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_keepdup-{keepdup}_no_control/{filler}/{sample}_peaks.broadPeak",
        xls =  "out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_keepdup-{keepdup}_no_control/{filler}/{sample}_peaks.xls"
    params:
        broad="--broad",
        format="--format {format}",
        gsize="--gsize {gsize}",
        keepdup="--keep-dup {keepdup}",
        name="--name {sample}",
        outdir="--outdir out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_keepdup-{keepdup}_no_control/{filler}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+",
        sample="[a-zA-Z0-9-_]+",
        keepdup="all|auto"
    shell:
        "{input.macs} callpeak -t {input.bam} {params}"


rule macs2_callpeak_broad_format_gsize_no_control_legacy:
    """
    Created:
        2018-02-10 00:14:58
    Aim:
        Straightforward peak calling, but not control by input...
        For ATAC:
        macs2 callpeak --treatment my.bam -g hs --bdg -q 0.05 --keep-dup all -f BAMPE --broad
    Test:
        out/macs2/callpeak_broad_format-BAM_gsize-hs_no_control/ln/alias/experiments/hg38_H3K27ac_thymus/input_peaks.broadPeak
    """
    input:
        macs = "opt/miniconda/envs/macs2/bin/macs2",
        bam =  "out/{filler}/{sample}.bam",
        bai =  "out/{filler}/{sample}.bam.bai"
    output:
        bed =  "out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_no_control/{filler}/{sample}_peaks.broadPeak",
        xls =  "out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_no_control/{filler}/{sample}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_no_control/{filler}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+",
        sample="[a-zA-Z0-9-_]+"
    shell:"""
    {input.macs} callpeak\
        --broad\
        -t {input.bam}\
        -f {wildcards.format}\
        -g {wildcards.gsize}\
        --name {wildcards.sample}\
        --outdir {params.outdir}
    """

rule macs2_callpeak_format_gsize_no_control_legacy:
    """
    Created:
        2018-02-10 00:13:48
    Aim:
        Straightforward peak calling, but not control by input...
    """
    input:
        macs = "opt/miniconda/envs/macs2/bin/macs2",
        bam =  "out/{filler}/{sample}.bam",
        bai =  "out/{filler}/{sample}.bam.bai"
    output:
        bed =  "out/macs2/callpeak_format-{format}_gsize-{gsize}_no_control/{filler}/{sample}_peaks.narrowPeak",
        xls =  "out/macs2/callpeak_format-{format}_gsize-{gsize}_no_control/{filler}/{sample}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_format-{format}_gsize-{gsize}_no_control/{filler}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+",
        sample="[a-zA-Z0-9-_]+"
    shell:"""
    {input.macs} callpeak\
        -t {input.bam}\
        -f {wildcards.format}\
        -g {wildcards.gsize}\
        --name {wildcards.sample}\
        --outdir {params.outdir}
    """

rule macs2_callpeak_no_control_format_gsize_call_summit_legacy:
    """
    """
    input:
        macs="opt/miniconda/envs/macs2/bin/macs2",
        bam="out/{sample}.bam",
        bai="out/{sample}.bam.bai"
    output:
        bed="out/macs2/callpeak_{format}_{gsize}_no_control_call_summit/{sample}.bed"
    params:
        outdir="out/macs2/callpeak_{format}_{gsize}_no_control_call_summit/{sample}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+"
    shell:"""
    {input.macs} callpeak -t {input.bam} \
        -f {wildcards.format} \
         -g {wildcards.gsize} \
        --call-summits \
        --outdir {params.outdir}

    touch {output.bed}
    """

#rule macs2_callpeak_broad_BAMPE:
#    """
#    """
#    input:
#        macs="opt/miniconda/envs/macs2/bin/macs2",
#        bam_chip="out/{id}/{chip}.bam",
#        bam_input="out/{id}/{input}.bam"
#    output:
#        bed="out/macs2/callpeak_broad_pe/{id}/{chip}_over_{input}_peaks.broadPeak"
#    params:
#        outdir="out/macs2/callpeak_broad_pe/{id}"
#    shell:"""
#    {input.macs} callpeak \
#        --treatment {input.bam_chip} \
#        --control {input.bam_input} \
#        --format BAMPE -g mm \
#        --name {wildcards.chip}_over_{wildcards.input} \
#        --outdir {params.outdir} \
#        --broad
#    """

rule macs2_callpeak_format_gsize_legacy:
    """
    Created:
        2017-07-17 11:34:32 - Use the new pattern for output files. I should update the other rules as well.
    Test:

    """
    input:
        macs="opt/miniconda/envs/macs2/bin/macs2",
        bam_chip="out/{id}/{chip}.bam",
        bam_input="out/{id}/{input}.bam"
    output:
        bed="out/macs2/callpeak_format-{format}_gsize-{gsize}/{id}/{chip}_over_{input}_peaks.narrowPeak",
        xls="out/macs2/callpeak_format-{format}_gsize-{gsize}/{id}/{chip}_over_{input}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_format-{format}_gsize-{gsize}/{id}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+"
    shell:
        """
        {input.macs} callpeak \
            --treatment {input.bam_chip} \
            --control {input.bam_input} \
            --format {wildcards.format} \
            --gsize {wildcards.gsize} \
            --name {wildcards.chip}_over_{wildcards.input} \
            --outdir {params.outdir}
        """

rule macs2_callpeak_format_gsize_nomodel_legacy:
    """
    Created:
        2018-03-16 19:11:11
    Aim:
        Blueprint settings for H3K27ac, H3K4me3, H3K9/14ac, H2A.Zac
    Note:
        macs2 callpeak -t chip.bam -n a_sensible_name --gsize hs -c input.bam --nomodel --shiftsize=half_fragment_size
    Test:
        out/macs2/callpeak_format-BAM_gsize-hs_nomodel_shiftsize-half_fragment_size/ln/alias/experiments/hg38_H3K27ac_thymus/CD34_over_input_peaks.narrowPeak
    """
    input:
        macs="opt/miniconda/envs/macs2/bin/macs2",
        bam_chip="out/{id}/{chip}.bam",
        bam_input="out/{id}/{input}.bam"
    output:
        bed    = "out/macs2/callpeak_format-{format}_gsize-{gsize}_nomodel/{id}/{chip}_over_{input}_peaks.narrowPeak",
        xls    = "out/macs2/callpeak_format-{format}_gsize-{gsize}_nomodel/{id}/{chip}_over_{input}_peaks.xls"
    params:
        outdir = "out/macs2/callpeak_format-{format}_gsize-{gsize}_nomodel/{id}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+",
        chip="[a-zA-Z0-9-_]+",
        input="[a-zA-Z0-9-_]+"
    shell:
        """
        {input.macs} callpeak \
            --treatment {input.bam_chip} \
            --control {input.bam_input} \
            --format {wildcards.format} \
            --gsize {wildcards.gsize} \
            --name {wildcards.chip}_over_{wildcards.input} \
            --nomodel\
            --outdir {params.outdir}
        """

rule macs2_callpeak_broad_format_gsize_legacy:
    """
    Created:
        2017-08-02 17:17:08
    Test:
        out/macs2/callpeak_broad_format-BAM_gsize-hs/ln/alias/experiments/hg38_H3K27ac_thymus/Th101_LC_over_input_peaks.broadPeak

    """
    input:
        macs="opt/miniconda/envs/macs2/bin/macs2",
        bam_chip="out/{id}/{chip}.bam",
        bam_input="out/{id}/{input}.bam"
    output:
        bed="out/macs2/callpeak_broad_format-{format}_gsize-{gsize}/{id}/{chip}_over_{input}_peaks.broadPeak",
        xls="out/macs2/callpeak_broad_format-{format}_gsize-{gsize}/{id}/{chip}_over_{input}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_broad_format-{format}_gsize-{gsize}/{id}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+",
        chip="[a-zA-Z0-9-_]+",
        input="[a-zA-Z0-9-_]+"
    shell:
        """
        {input.macs} callpeak \
            --broad \
            --treatment {input.bam_chip} \
            --control {input.bam_input} \
            --format {wildcards.format} \
            --gsize {wildcards.gsize} \
            --name {wildcards.chip}_over_{wildcards.input} \
            --outdir {params.outdir}
        """

rule macs2_callpeak_format_gsize_SPMR_legacy:
    """
    Created:
        2017-07-17 11:34:32 - Use the new pattern for output files. I should update the other rules as well.
    Test:

    """
    input:
        macs="opt/miniconda/envs/macs2/bin/macs2",
        bam_chip="out/{id}/{chip}.bam",
        bam_input="out/{id}/{input}.bam"
    output:
        bed="out/macs2/callpeak_format-{format}_gsize-{gsize}_SPMR/{id}/{chip}_over_{input}_peaks.narrowPeak",
        xls="out/macs2/callpeak_format-{format}_gsize-{gsize}_SPMR/{id}/{chip}_over_{input}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_format-{format}_gsize-{gsize}_SPMR/{id}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+"
    shell:
        """
        {input.macs} callpeak \
            --treatment {input.bam_chip} \
            --control {input.bam_input} \
            --format {wildcards.format} \
            --gsize {wildcards.gsize} \
            --name {wildcards.chip}_over_{wildcards.input} \
            --outdir {params.outdir} \
            --SPMR
        """

rule macs2_callpeak_broad_format_gsize_SPMR_legacy:
    """
    Created:
        2017-10-05 09:55:16
    Test:

    """
    input:
        macs=     "opt/miniconda/envs/macs2/bin/macs2",
        bam_chip= "out/{id}/{chip}.bam",
        bam_input="out/{id}/{input}.bam"
    output:
        bed=   "out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_SPMR/{id}/{chip}_over_{input}_peaks.broadPeak",
        xls=   "out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_SPMR/{id}/{chip}_over_{input}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_broad_format-{format}_gsize-{gsize}_SPMR/{id}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+"
    shell:
        """
        {input.macs} callpeak \
            --treatment {input.bam_chip} \
            --control {input.bam_input} \
            --broad \
            --format {wildcards.format} \
            --gsize {wildcards.gsize} \
            --name {wildcards.chip}_over_{wildcards.input} \
            --outdir {params.outdir} \
            --SPMR
        """

#rule macs2_callpeak_broad_format_gsize_SPMR_legacy:
#    """
#    Modified:
#        2017-07-13 12:04:29 - Changed to legacy because of the way SPMR argument is used with 'True' and 'False'.
#    This argument may be valuable to compute FRIP-like metric:
#    --SPMR  If True, MACS will save signal per million reads for
#            fragment pileup profiles. Require -B to be set.
#            Default: False
#    https://github.com/taoliu/MACS/issues/55#issuecomment-58199145
#
#    Test:
#        out/macs2/callpeak_broad_BAMPE_mm_SPMRTrue/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.xls
#    """
#    input:
#        macs="opt/miniconda/envs/macs2/bin/macs2",
#        bam_chip="out/{id}/{chip}.bam",
#        bam_input="out/{id}/{input}.bam"
#    output:
#        bed="out/macs2/callpeak_broad_{format}_{gsize}_SPMR{spmr}/{id}/{chip}_over_{input}_peaks.broadPeak",
#        xls="out/macs2/callpeak_broad_{format}_{gsize}_SPMR{spmr}/{id}/{chip}_over_{input}_peaks.xls"
#    params:
#        outdir="out/macs2/callpeak_broad_{format}_{gsize}_SPMR{spmr}/{id}"
#    wildcard_constraints:
#        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
#        gsize="mm|hs|ce|dm|[0-9]+",
#        spmr="True|False"
#    shell:"""
#    if [ '{wildcards.spmr}' == 'True' ];
#    then
#        {input.macs} callpeak \
#            --treatment {input.bam_chip} \
#            --control {input.bam_input} \
#            --format {wildcards.format} --gsize {wildcards.gsize} \
#            --name {wildcards.chip}_over_{wildcards.input} \
#            --outdir {params.outdir} \
#            --broad \
#            --bdg \
#            --SPMR
#    else
#        {input.macs} callpeak \
#            --treatment {input.bam_chip} \
#            --control {input.bam_input} \
#            --format {wildcards.format} --gsize {wildcards.gsize} \
#            --name {wildcards.chip}_over_{wildcards.input} \
#            --outdir {params.outdir} \
#            --broad
#    fi
#    """

rule macs2_callpeak_broad_qvalue_format_gsize_legacy:
    """
    Created: 2017
    MACS peak calling.

    -q QVALUE, --qvalue QVALUE
                        Minimum FDR (q-value) cutoff for peak detection.
                        DEFAULT: 0.05. -q, and -p are mutually exclusive
    Test:
        out/macs2/callpeak_broad_qvalue0.01_BAMPE_mm/bam/mm9/merge_run140_run141/H4K5ac-Nut-WT_over_Input-Nut-WT_peaks.xls
    """
    input:
        macs="opt/miniconda/envs/macs2/bin/macs2",
        bam_chip="out/{id}/{chip}.bam",
        bam_input="out/{id}/{input}.bam"
    output:
        bed="out/macs2/callpeak_broad_qvalue{qvalue}_{format}_{gsize}/{id}/{chip}_over_{input}_peaks.broadPeak",
        xls="out/macs2/callpeak_broad_qvalue{qvalue}_{format}_{gsize}/{id}/{chip}_over_{input}_peaks.xls"
    params:
        outdir="out/macs2/callpeak_broad_qvalue{qvalue}_{format}_{gsize}/{id}"
    wildcard_constraints:
        format="BAM|BAMPE|AUTO|SAM|BED|ELAND|ELANDMULTI|ELANDEXPORT|BOWTIE|BEDPE",
        gsize="mm|hs|ce|dm|[0-9]+",
        qvalue="0\.[0-9]+"
    shell:"""
    {input.macs} callpeak \
        --treatment {input.bam_chip} \
        --control {input.bam_input} \
        --format {wildcards.format} --gsize {wildcards.gsize} \
        --name {wildcards.chip}_over_{wildcards.input} \
        --outdir {params.outdir} \
        --broad \
        --qvalue {wildcards.qvalue}
    """
