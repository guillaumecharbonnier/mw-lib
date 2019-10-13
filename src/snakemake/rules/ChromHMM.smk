def input_bed_ChromHMM_BinarizeBed_dependencies(wildcards):
    """
    Created:
        2017-06-14 16:42:12
    Aim:
        A function retrieving the third column of 'cellmarkfiletable' to return the list of bed files needed by ChromHMM BinarizeBed and ensure bed files are created by Snakemake before being used by ChromHMM.
    """
    cellmarkfiletable = config['ids'][wildcards['cellmarkfiletable_id']]

    df = pandas.read_table(cellmarkfiletable, header=None, comment='#')
    # Column 3 should contain samples
    paths=df[2].tolist()
    # Column 4 may contain inputs
    if len(df.columns) ==4:
        paths = paths + df[3].tolist()

    return(paths)

rule ChromHMM_BinarizeBed_cellmarkfiletable_chrominfo_extra:
    """
    Created:
        2018-05-02 23:18:38
    Aim:
        General purpose binarisation function.
        #IDée: archive le contenu de l'output, déclarer l'archive, puis open l'archive dans la règle suivante.
    Doc:
        -b binsize  The number of base pairs in a bin determining the resolution of the model learning and segmentation.  By default this parameter value is set to 200base pairs
    Test:
        out/ChromHMM/BinarizeBed_b-200_chrominfo-hg19-main-chr_merged-Blueprint-thymic-populations-with-input-as-control-hg19/done
        out/ChromHMM/BinarizeBed_b-200_chrominfo-hg38-main-chr_Blueprint-thymic-populations-with-input-as-control/done
        out/ChromHMM/developBinarizeExtra_chrominfo-mm10_cellmarkfiletable-test-srr.tar.gz
    Doc:
        usage BinarizeBed [-b binsize][-c controldir][-center][-colfields chromosome,start,end[,strand]][-e offsetend][-f foldthresh][-g signalthresh][-n shift][-o outputcontroldir][-p poissonthresh][-peaks][-s offsetstart][-strictthresh][-t outputsignaldir][-u pseudocountcontrol][-w flankwidthcontrol] chromosomelengthfile inputbeddir cellmarkfiletable outputbinarydir
    """
    input:
        chromosomelengthfile = lambda wildcards: config['ids'][wildcards.chrominfo_id],
        cellmarkfiletable = lambda wildcards: config['ids'][wildcards.cellmarkfiletable_id],
        bed = input_bed_ChromHMM_BinarizeBed_dependencies
    output:
        done=touch("out/{tool}{extra}_{chrominfo_id}_{cellmarkfiletable_id}/done")
        #binarizedbedtar = "out/{tool}{extra}_{chrominfo_id}_{cellmarkfiletable_id}.tar.gz"
    log:
        "out/{tool}{extra}_{chrominfo_id}_{cellmarkfiletable_id}/log"
    benchmark:
        "out/{tool}{extra}_{chrominfo_id}_{cellmarkfiletable_id}/benchmark.tsv"
    params:
        inputbeddir=".",
        outputbinarydir="out/{tool}{extra}_{chrominfo_id}_{cellmarkfiletable_id}",
        extra = params_extra
    wildcard_constraints:
        tool="ChromHMM/BinarizeBed",
        cellmarkfiletable_id="[\w-]+",
    conda:
        "../envs/chromhmm.yaml"
    shell:
        """
        ChromHMM.sh\
            -Xmx32768m\
            BinarizeBed\
            {params.extra}\
            {input.chromosomelengthfile}\
            {params.inputbeddir}\
            {input.cellmarkfiletable}\
            {params.outputbinarydir} &> {log}
        """
        #tar zcvf {output.binarizedbedtar} -C {params.outputbinarydir} .


#rule test_aggregate_checkpoint:
#    input:
#        checkpoints.ChromHMM_BinarizeBed_cellmarkfiletable_chrominfo_extra.get(tool="ChromHMM/developBinarizeExtra", extra=
#    output:
#        "test_aggregate_checkpoint.txt"
#    shell:
#        "touch {output}"
#        

rule ChromHMM_LearnModel_extra:
    """
    Created:
        2019-09-04 17:47:48
    Test:
        out/ChromHMM/LearnModel_numstates-2_assembly-mm10/ChromHMM/developBinarizeExtra_chrominfo-mm10_cellmarkfiletable-test-srr/model_2.txt
    """
    input:
        binarizedbedtar = "out/{filler}.tar.gz"
    output:
        #done = touch("out/{tool}{extra}_numstates-{numstates}_assembly-{assembly}/{filler}/done"),
        trans_emiss = expand("out/{{tool}}{{extra}}_numstates-{{numstates}}_assembly-{{assembly}}/{{filler}}/{trans_emiss}_{{numstates}}.{ext}", trans_emiss = ["transitions","emissions"], ext=["png","svg","txt"]),
        model       = "out/{tool}{extra}_numstates-{numstates}_assembly-{assembly}/{filler}/model_{numstates}.txt",
        #webpage     = "out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/webpage_{numstates}.html",
        # If any stage-dependent file are needed for downstream analysis, they have to be explicitely defined in "ln_gather_dynamic" rule
    params:
        extra = params_extra
    wildcard_constraints:
        numstates="[0-9]+",
        assembly="mm10|hg38",
        tool="ChromHMM/LearnModel" #|_custom_features" Add custom features if needed 
#    benchmark:
#        "log/snakemake/benchmark/out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}.txt"
    conda:
        "../envs/chromhmm.yaml"
    threads:
        16
    shell:
        """
        OUTPUTDIR=`dirname {output.model}`
        rm -rf $OUTPUTDIR/input
        mkdir $OUTPUTDIR/input
        tar zxvf {input.binarizedbedtar} -C $OUTPUTDIR/input

        ChromHMM.sh -Xmx32768m \
            LearnModel \
            -p {threads} \
            $OUTPUTDIR/input \
            $OUTPUTDIR \
            {wildcards.numstates} \
            {wildcards.assembly}
        """

rule ChromHMM_MakeBrowserFiles:
    """
    Created:
        2018-07-26 11:56:21
    Aim:
        Create files to browse states in IGV
    Test:
    sm\
        out/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/CD34_11_segments.done\
        out/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/EC_11_segments.done\
        out/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/LC_11_segments.done\
        out/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/SP4_11_segments.done\
        out/ChromHMM/MakeBrowserFiles/ChromHMM/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/SP8_11_segments.done
    """
    input:
        segments = "out/{filler}/{filename}.bed"
    output:
        dense    = "out/ChromHMM/MakeBrowserFiles/{filler}/{filename}_dense.bed",
        expanded = "out/ChromHMM/MakeBrowserFiles/{filler}/{filename}_expanded.bed"
    params:
        memory="65536m" # previously 32768 but not enough for test9
    conda:
        "../envs/chromhmm.yaml"
    shell:
        """
        outputfileprefix="out/ChromHMM/MakeBrowserFiles/{wildcards.filler}/{wildcards.filename}"

        ChromHMM.sh\
            -Xmx{params.memory}\
            MakeBrowserFiles\
            {input.segments}\
            {wildcards.filename}\
            $outputfileprefix
        """


#rule ChromHMM_MakeSegments_model_binarization_done:
#    """
#    Created:
#        2018-05-10 20:13:59
#    Note:
#        out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/
#    """
#    input:
#        model="out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/model_11_All_cell_types.txt",
#        binarization_done="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}/done"
#    output:
#        done="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/done",
#        png_emission="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/emissions_{numstates}.png",
#        png_transition="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/transitions_{numstates}.png",
#        #bed_segments="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/PARSETHETABLEFORSTAGEHERE_{numstates}.png"
#    params:
#        inputdir="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}",
#        outputdir="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}"
#    params:
#        inpdir="out/ChromHMM/BinarizeBed_b-200_chrominfo-hg19-main-chr_Blueprint-thymic-populations-with-input-as-control-hg19",
#        outdir="out/ChromHMM/MakeSegments_our_thymic_samples_into_chromdet_original_paper",
#        memory="65536m" # previously 32768 but not enough for test9
#    shell:
#        """
#        {input.chromhmm}\
#            -Xmx{params.memory}\
#            MakeSegmentation\
#            {input.model}\
#            {params.inpdir}\
#            {params.outdir}
#        """
#



# Legacy rules below
#
#rule ChromHMM_BinarizeBed_b_chrominfo_broken:
#    """
#    BROKEN by_the_rewritten_input_bed_ChromHMM_BinarizeBed_dependencies_function:
#    
#    Created:
#        2018-05-02 23:18:38
#    Aim:
#        General purpose binarisation function.
#    Test:
#        out/ChromHMM/BinarizeBed_b-200_chrominfo-hg19-main-chr_merged-Blueprint-thymic-populations-with-input-as-control-hg19/done
#        out/ChromHMM/BinarizeBed_b-200_chrominfo-hg38-main-chr_Blueprint-thymic-populations-with-input-as-control/done
#    """
#    input:
#        #chromhmm="opt/miniconda/envs/chromhmm/bin/ChromHMM.sh",
#        chromosomelengthfile = lambda wildcards: config['ids'][wildcards.chrominfo_id],
#        cellmarkfiletable="src/chromhmm/cellmarkfiletable/{cellmarkfiletable_id}.tsv",
#        bed=input_bed_ChromHMM_BinarizeBed_dependencies
#    output:
#        done=touch("out/ChromHMM/BinarizeBed_b-{binSize}_{chrominfo_id}_{cellmarkfiletable_id}/done")
#    params:
#        inputbeddir=".",
#        outputbinarydir="out/ChromHMM/BinarizeBed_b-{binSize}_chrominfo-{chrominfo_id}_{cellmarkfiletable_id}"
#    wildcard_constraints:
#        cellmarkfiletable_id="[a-zA-Z0-9-]+",
#        binSize="[0-9]+"
#    conda:
#        "../envs/chromhmm.yaml"
#    shell:
#        """
#        ChromHMM.sh\
#            -Xmx32768m\
#            BinarizeBed\
#            -b {wildcards.binSize}\
#            -center\
#            {input.chromosomelengthfile}\
#            {params.inputbeddir}\
#            {input.cellmarkfiletable}\
#            {params.outputbinarydir}
#        """

rule ChromHMM_BinarizeBed:
    """
    Created:
        2017-06-14 11:20:22
    Aim:
        Trying using wdir for inputbeddir.
        This can allow to use a function to parse 'cellmarkfiletable' for the paths of bed files needed.
    Test:
        "out/ChromHMM/BinarizeBed_test2/done"
    """
    input:
        chromosomelengthfile="out/awk/extract_main_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt",
        cellmarkfiletable="src/chromhmm/cellmarkfiletable/{cellmarkfiletable_id}.tsv",
        bed = input_bed_ChromHMM_BinarizeBed_dependencies
    output:
        done="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}/done"
    params:
        inputbeddir=".",
        outputbinarydir="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}"
    wildcard_constraints:
        cellmarkfiletable_id="[a-zA-Z0-9-]+"
    conda:
        "../envs/chromhmm.yaml"
    shell:
        """
        ChromHMM.sh -Xmx32768m \
            BinarizeBed \
            -b 200 \
            -center \
            {input.chromosomelengthfile} \
            {params.inputbeddir} \
            {input.cellmarkfiletable} \
            {params.outputbinarydir}

        touch {output.done}
        """

rule ChromHMM_BinarizeBed_b:
    """
    Created:
        2017-11-17 19:22:46
    Aim:
        I want to try to use bin size smaller than 200~bp.
    Test:
        "out/ChromHMM/BinarizeBed_b-50_test19/done"
    """
    input:
        chromosomelengthfile="out/awk/extract_main_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt",
        cellmarkfiletable="src/chromhmm/cellmarkfiletable/{cellmarkfiletable_id}.tsv",
        bed=input_bed_ChromHMM_BinarizeBed_dependencies
    output:
        done="out/ChromHMM/BinarizeBed_b-{binSize}_{cellmarkfiletable_id}/done"
    params:
        inputbeddir=".",
        outputbinarydir="out/ChromHMM/BinarizeBed_b-{binSize}_{cellmarkfiletable_id}"
    wildcard_constraints:
        cellmarkfiletable_id="[a-zA-Z0-9-]+",
        binSize="[0-9]+"
    conda:
        "../envs/chromhmm.yaml"
    shell:
        """
        ChromHMM.sh -Xmx32768m\
            BinarizeBed\
            -b {wildcards.binSize}\
            -center\
            {input.chromosomelengthfile}\
            {params.inputbeddir}\
            {input.cellmarkfiletable}\
            {params.outputbinarydir}

        touch {output.done}
        """

rule ChromHMM_LearnModel_numstates_assembly:
    """
    Created:
        2017-06-13 10:47:36
    Test:
        out/ChromHMM/LearnModel_test1_numstates-10_assembly-mm10/done
        out/ChromHMM/LearnModel_test3_numstates-10_assembly-mm10/done
        out/ChromHMM/LearnModel_test4_numstates-20_assembly-mm10/done
        out/ChromHMM/LearnModel_test13_numstates-20_assembly-mm10/done
        out/ChromHMM/LearnModel_test19_numstates-10_assembly-mm10/done
    """
    input:
        chromhmm="opt/miniconda/envs/chromhmm{_custom_features}/bin/ChromHMM.sh",
        chromosomelengthfile="out/awk/extract_main_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt",
        done="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}/done"
    output:
        done = touch("out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/done"),
        trans_emiss = expand("out/ChromHMM{{_custom_features}}/LearnModel_{{cellmarkfiletable_id}}_numstates-{{numstates}}_assembly-{{assembly}}/{trans_emiss}_{{numstates}}.{ext}", trans_emiss = ["transitions","emissions"], ext=["png","svg","txt"]),
        model       = "out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/model_{numstates}.txt",
        webpage     = "out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/webpage_{numstates}.html",
        # If any stage-dependent file are needed for downstream analysis, they have to be explicitely defined in "ln_gather_dynamic" rule
    params:
        inputdir="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}",
        outputdir="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}"
    wildcard_constraints:
        numstates="[0-9]+",
        assembly="mm10|hg38",
        _custom_features="|_custom_features"
    benchmark:
        "log/snakemake/benchmark/out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}.txt"
    threads:
        16
    shell:
        """
        {input.chromhmm} -Xmx32768m \
            LearnModel \
            -p {threads} \
            {params.inputdir} \
            {params.outputdir} \
            {wildcards.numstates} \
            {wildcards.assembly}
        """

rule ChromHMM_LearnModel_numstates_assembly_seed:
    """
    Created:
        2017-06-20 13:49:58
    Aim:
        Seed is added to assess the reproducibility of the learning.
    Test:
        out/ChromHMM/LearnModel_test7_numstates-40_assembly-mm10_seed-1/done
    """
    input:
        chromhmm="opt/miniconda/envs/chromhmm/bin/ChromHMM.sh",
        chromosomelengthfile="out/awk/extract_main_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt",
        done="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}/done"
    output:
        done="out/ChromHMM/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}_seed-{seed}/done",
        png_emission="out/ChromHMM/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}_seed-{seed}/emissions_{numstates}.png",
        png_transition="out/ChromHMM/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}_seed-{seed}/transitions_{numstates}.png",
        #bed_segments="out/ChromHMM/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/PARSETHETABLEFORSTAGEHERE_{numstates}.png"
    params:
        inputdir="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}",
        outputdir="out/ChromHMM/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}_seed-{seed}",
        memory="65536m" # previously 32768 but not enough for test9
    wildcard_constraints:
        numstates="[0-9]+",
        assembly="mm10|hg38",
        seed="[0-9]+"
    benchmark:
        "log/snakemake/benchmark/out/ChromHMM/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}_seed-{seed}.txt"
    threads:
        16
    shell:
        """
        {input.chromhmm} -Xmx{params.memory} \
            LearnModel \
            -p {threads} \
            -s {wildcards.seed} \
            {params.inputdir} \
            {params.outputdir} \
            {wildcards.numstates} \
            {wildcards.assembly}

        touch {output.done}
        """


rule ChromHMM_MakeSegments_model_binarization_done:
    """
    Created:
        2018-05-10 20:13:59
    Note:
        out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/
    Test:
        out/MakeSegments_our_merged_thymic_samples_into_chromdet_original_paper/EC_11_segments_dense.bed
    """
    input:
        model="out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/model_11_All_cell_types.txt",
        binarization_done="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}/done"
    output:
        done="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/done",
        png_emission="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/emissions_{numstates}.png",
        png_transition="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/transitions_{numstates}.png",
        #bed_segments="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}/PARSETHETABLEFORSTAGEHERE_{numstates}.png"
    params:
        inputdir="out/ChromHMM/BinarizeBed_{cellmarkfiletable_id}",
        outputdir="out/ChromHMM{_custom_features}/LearnModel_{cellmarkfiletable_id}_numstates-{numstates}_assembly-{assembly}"
    params:
        inpdir="out/ChromHMM/BinarizeBed_b-200_chrominfo-hg19-main-chr_Blueprint-thymic-populations-with-input-as-control-hg19",
        outdir="out/ChromHMM/MakeSegments_our_thymic_samples_into_chromdet_original_paper",
        memory="65536m" # previously 32768 but not enough for test9
    conda:
        "../envs/chromhmm.yaml"
    shell:
        """
        ChromHMM.sh\
            -Xmx{params.memory}\
            MakeSegmentation\
            {input.model}\
            {params.inpdir}\
            {params.outdir}
        """


rule ChromHMM_BinarizeBed_test1:
    """
    Created:
        2017-06-13 10:47:36
    """
    input:
        chromhmm="opt/miniconda/envs/chromhmm/bin/ChromHMM.sh",
        chromosomelengthfile="out/awk/extract_main_chr/gunzip/to-stdout/rsync/ucsc/goldenPath/mm10/database/chromInfo.txt",
        cellmarkfiletable="src/chromhmm/cellmarkfiletable/test1.tsv"
    output:
    params:
        inputbeddir="out/ln/make_ChromHMM_inputbeddir_test1",
        outputbinarydir="out/ChromHMM/BinarizeBed_test1"
    shell:
        """
        {input.chromhmm} \
            BinarizeBed \
            -b 200 \
            -center \
            {input.chromosomelengthfile} \
            {params.inputbeddir} \
            {input.cellmarkfiletable} \
            {params.outputbinarydir}
        """

rule ln_make_ChromHMM_input_test1:
    input:
        band1="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.bed",
        band2="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.bed",
        band3="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band3.bed",
        band4="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-100_lmax-130/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band4.bed",
        band5="out/awk/extract_main_chr/awk/convert_bedpe_to_bed6_insert_size/bedtools/bamtobed_bedpe/samtools/sort_n/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band5.bed"
    output:
        bed=expand("out/ln/make_ChromHMM_inputbeddir_test1/MNase_Spm_WT_band{band}.bed", band=["1","2","3","4","5"])
    params:
        inputbeddir="out/ln/make_ChromHMM_inputbeddir_test1"
    shell:
        """
        ln {input.band1} {params.inputbeddir}/MNase_Spm_WT_band1.bed
        ln {input.band2} {params.inputbeddir}/MNase_Spm_WT_band2.bed
        ln {input.band3} {params.inputbeddir}/MNase_Spm_WT_band3.bed
        ln {input.band4} {params.inputbeddir}/MNase_Spm_WT_band4.bed
        ln {input.band5} {params.inputbeddir}/MNase_Spm_WT_band5.bed
        """

rule tmp_test1_bw:
    """
    Created:
        2017-06-13 15:56:06
    Aim:
        Check the states called by ChromHMM correspond to something we can see with coverage tracks.
    """
    input:
        band1="out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.bw",
        band2="out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.bw",
        band3="out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band3.bw",
        band4="out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-100_lmax-130/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band4.bw",
        band5="out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band5.bw",
        band6="out/deepTools/2_5_1_bamCoverage_binSize-10_minMappingQuality-0_normalizeUsingRPKM_extendReads-/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band6.bw"


rule tmp_test19_bw:
    """
    Created:
        2017-11-17 13:08:56
    Aim:
        Check the states called by ChromHMM correspond to something we can see with coverage tracks.
    """
    input:
        spermatozoa_mnase_nuc="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184_run187/s1-MNase_Spm_WT_band1_s2-MNase_Spm_WT_band3_s3-MNase_Spm_WT_band6.bw",
        spermatozoa_mnase_ss="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge_three_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167_run184/s1-MNase_Spm_WT_band2_s2-MNase_Spm_WT_band4_s3-MNase_Spm_WT_band5.bw",
        pachytene_mnase_nuc="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/s1-MNS-P-WT_s2-MNS-P-KO.bw",
        round_spermatid_mnase_nuc="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run119_run124/s1-MNS-R-WT_s2-MNS-R-KO.bw",
        condensing_spermatid_mnase_nuc="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/s1-MNS-SC-WT_s2-MNS-SC-KO.bw",
        condensing_spermatid_mnase_ss="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/merge_two_samples/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run113_run125_run126/s1-PSK-SC-WT_s2-PSK-SC-KO.bw",
        naked_dna_nuc_like="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/merge_three_samples/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/ln/rename_run192_tgml/s1-mnase_naked_dna_1_s2-mnase_naked_dna_2_s3-mnase_naked_dna_3.bw",
        naked_dna_ss_like="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/merge_three_samples/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/ln/rename_run192_tgml/s1-mnase_naked_dna_4_s2-mnase_naked_dna_5_s3-mnase_naked_dna_6.bw",
        naked_dna_true_nuc_like="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge_three_samples/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/ln/rename_run192_tgml/s1-mnase_naked_dna_1_s2-mnase_naked_dna_2_s3-mnase_naked_dna_3.bw",
        naked_dna_true_ss_like="out/deepTools/bamCoverage_binSize-10_normalizeUsingRPKM_extendReads-/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge_three_samples/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/ln/rename_run192_tgml/s1-mnase_naked_dna_4_s2-mnase_naked_dna_5_s3-mnase_naked_dna_6.bw"


rule tmp_test1_idxstats:
    """
    Created:
        2017-06-13 15:56:06
    Aim:
        Check the states called by ChromHMM correspond to something we can see with coverage tracks.
    """
    input:
        band1="out/samtools/idxstats/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.tsv",
        band2="out/samtools/idxstats/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.tsv",
        band3="out/samtools/idxstats/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band3.tsv",
        band4="out/samtools/idxstats/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-100_lmax-130/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band4.tsv",
        band5="out/samtools/idxstats/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band5.tsv"

rule tmp_test1_flagstat:
    """
    Created:
        2017-06-13 17:13:03
    Aim:
        Check the states called by ChromHMM correspond to something we can see with coverage tracks.
    """
    input:
        band1="out/samtools/flagstat/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band1.tsv",
        band2="out/samtools/flagstat/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band2.tsv",
        band3="out/samtools/flagstat/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-130_lmax-170/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band3.tsv",
        band4="out/samtools/flagstat/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-100_lmax-130/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band4.tsv",
        band5="out/samtools/flagstat/samtools/sort/samtools/view_bSh/java/select_subpopulations_from_bam_lmin-30_lmax-100/samtools/merge/samtools/sort/samtools/view_bSh/bowtie2/pe_mm10/ln/paired_end_remove_mate_prefix/fastx_toolkit/fastx_trimmer_l-30/gunzip/merge_lanes_nextseq500_paired_end/inp/fastq/run163_run167/MNase_Spm_WT_band5.tsv"
