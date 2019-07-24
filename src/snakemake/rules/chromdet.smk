rule chromdet_alternative_collapse_segments:
    """
    Created:
        2018-06-18 12:02:31
    Aim:
        I have realized States_collapse from github merge "transcription" and "enhancer"
    Test:
        out/chromdet/alternative-collapse_segments-hg19-segments-our-complete-merged-thymic-samples-into-chromdet-original-paper-t-cells/samples_collapsed_filtered_chromatin_space.tsv out/chromdet/alternative-collapse_segments-hg19-segments-our-complete-merged-thymic-samples-into-chromdet-original-paper/samples_collapsed_filtered_chromatin_space.tsv
    """
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        segments = lambda wildcards: eval(config['ids'][wildcards.bed_list_id]),
        annotation="src/chromdet/annotation.tsv",
        #states_collapse="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master/ChromDet-master/test/States_collapse.txt"
        states_collapse="src/chromdet/states_collapse.tsv"
    output:
        expand("out/chromdet/alternative-collapse_segments-{{bed_list_id}}/samples_collapsed{file}",file=CHROMDET_OUTPUT_FILES)
    params:
        outdir="out/chromdet/alternative-collapse_segments-{bed_list_id}"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        #export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        #export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        #export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        # This is ugly but has to be done because some commands in chromdet are relative paths.
        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -a {WDIR}/{input.annotation} -c {WDIR}/{input.states_collapse} -s {WDIR}/opt/miniconda/envs/chromdet/bin/S3Det_modified/ -v
        
        #cd {input.dir}/ChromDet-master/scripts
        #run_S3det_analysis.pl -d {params.outdir} -a {input.annotation} -c {input.states_collapse} -s {WDIR}/opt/miniconda/envs/chromdet/bin/S3Det_modified/ -v

        """

rule chromdet_segments:
    """
    Created:
        2018-05-03 16:08:23
    Aim:
        General purpose chromdet rule.
    Test:
        out/chromdet/segments-hg19-segments-our-complete-merged-thymic-samples-into-chromdet-original-paper-t-cells/samples_collapsed_filtered_chromatin_space.tsv
        out/chromdet/segments-hg19-segments-our-complete-merged-thymic-samples-into-chromdet-original-paper/samples_collapsed_filtered_chromatin_space.tsv
        out/chromdet/segments-hg19-segments-our-complete-thymic-samples-into-chromdet-original-paper/samples_collapsed_filtered_chromatin_space.tsv
        out/chromdet/segments-hg19-segments-our-thymic-samples-into-chromdet-original-paper/samples_collapsed_filtered_chromatin_space.tsv
    """
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        segments = lambda wildcards: eval(config['ids'][wildcards.bed_list_id]),
        annotation="src/chromdet/annotation.tsv",
        states_collapse="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master/ChromDet-master/test/States_collapse.txt"
    output:
        expand("out/chromdet/segments-{{bed_list_id}}/samples_collapsed{file}",file=CHROMDET_OUTPUT_FILES)
    params:
        outdir="out/chromdet/segments-{bed_list_id}"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        #export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        #export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        #export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        # This is ugly but has to be done because some commands in chromdet are relative paths.
        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -a {WDIR}/{input.annotation} -c {WDIR}/{input.states_collapse} -s {WDIR}/opt/miniconda/envs/chromdet/bin/S3Det_modified/ -v
        
        #cd {input.dir}/ChromDet-master/scripts
        #run_S3det_analysis.pl -d {params.outdir} -a {input.annotation} -c {input.states_collapse} -s {WDIR}/opt/miniconda/envs/chromdet/bin/S3Det_modified/ -v

        """

rule chromdet_no_collapse_segments:
    """
    Created:
        2018-07-28 00:39:35
    Aim:
        I need to get the matrix from samples_collapsed 
    Test:
        out/chromdet/no-collapse_segments-hg19-segments-our-complete-merged-thymic-samples-into-chromdet-original-paper-mhsc-t-cells-filt/samples_collapsed_filtered_chromatin_space.tsv out/chromdet/no-collapse_segments-hg19-segments-our-complete-merged-thymic-samples-into-chromdet-original-paper/samples_collapsed_filtered_chromatin_space.tsv
    """
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        segments = lambda wildcards: eval(config['ids'][wildcards.bed_list_id]),
        annotation="src/chromdet/annotation.tsv",
        #states_collapse="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master/ChromDet-master/test/States_collapse.txt"
        #states_collapse="src/chromdet/states_collapse.tsv"
    output:
        expand("out/chromdet/no-collapse_segments-{{bed_list_id}}/samples_collapsed{file}",file=CHROMDET_OUTPUT_FILES)
    params:
        outdir="out/chromdet/no-collapse_segments-{bed_list_id}"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u

        # This is ugly but has to be done because some commands in chromdet are relative paths.
        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -a {WDIR}/{input.annotation} -s {WDIR}/opt/miniconda/envs/chromdet/bin/S3Det_modified/ -v
        """


rule chromdet_toy_example:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master"
    shell:
        """
        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d ../test/ -c ../test/States_collapse.txt -a ../test/Samples_beds.tsv -s {WDIR}/opt/S3Det_modified/ -v
        """

rule chromdet_toy_example_conda:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        chromdet="opt/miniconda/envs/chromdet/bin/run_S3det_analysis.pl"
    params:
        S3Det_dir="opt/miniconda/envs/chromdet/bin/S3Det_modified/"
    shell:
        """
        set +u; source activate chromdet_after_patching_ln_regex; set -u

        run_S3det_analysis.pl -d ../test/ -c ../test/States_collapse.txt -a ../test/Samples_beds.tsv -s {WDIR}/{params.S3Det_dir} -v
        """


rule tar_xvzf_chromdet_redo_original_paper:
    """
    Created:
        2018-04-13 13:06:22
    Aim:
        Original paper has the amazing idea to work with hematopoietic cells.
        Let's reproduce its main MCA plot, then I will add our samples.
        Beware GRCh37 and not 38.
    """
    input:
        files=expand("out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/{file}",file=CHROMATIN_STATES_CARRILLO_BUILD37),
        archive=expand("out/wget/ftp/ftp.ebi.ac.uk/pub/databases/blueprint/paper_data_sets/chromatin_states_carrillo_build37/{archive}.tar.gz", archive=["BLUEPRINT_cell_lines","BLUEPRINT_disease","BLUEPRINT_healthy","ENCODE","ROADMAP"])
    output:
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER)
    params:
        outdir="out/tar/xvzf_chromdet_redo_original_paper"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for ARCHIVE in {input.archive}
        do
            echo "Extracting $ARCHIVE"
            tar -xvzf $ARCHIVE --directory {params.outdir}
        done
        """

rule ln_chromdet_redo_original_paper:
    input:
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER)
    output:
        segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER)
    params:
        outdir="out/ln/chromdet_redo_original_paper"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done
        """

rule chromdet_redo_original_paper:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER)
    output:
        "out/ln/chromdet_redo_original_paper/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        indir="out/ln/chromdet_redo_original_paper",
        outdir="out/chromdet/redo_original_paper"
    shell:
        """
        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.indir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v

        """

rule chromdet_redo_original_paper_blueprint_healthy:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER)
    output:
        "out/tar/xvzf_chromdet_redo_original_paper/BLUEPRINT_healthy/SEGMENTATION/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        indir="out/tar/xvzf_chromdet_redo_original_paper/BLUEPRINT_healthy/SEGMENTATION"
        #outdir="out/chromdet/redo_original_paper"
    shell:
        """
        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.indir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v

        """

rule chromdet_redo_original_paper_table_s1:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
    output:
        "out/chromdet/redo_original_paper_table_s1/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v

        """

rule chromdet_redo_original_paper_table_s1_without_collapse:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
    output:
        "out/chromdet/redo_original_paper_table_s1_without_collapse/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1_without_collapse"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -s {WDIR}/opt/S3Det_modified/ -v

        """

rule chromdet_redo_original_paper_table_s1_with_disease:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1_WITH_DISEASE)
    output:
        "out/chromdet/redo_original_paper_table_s1_with_disease/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1_with_disease"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v

        """


rule chromdet_redo_original_paper_table_s1_with_disease_without_collapse:
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1_WITH_DISEASE)
    output:
        "out/chromdet/redo_original_paper_table_s1_with_disease_without_collapse/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1_with_disease_without_collapse"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -s {WDIR}/opt/S3Det_modified/ -v

        """

rule chromdet_redo_original_paper_table_s1_autosomes:
    """
    Created:
        2018-04-30 11:32:10
    Aim:
        segmentation files were containing all chromosomes whereas material and methods explicitly use only autosomes. This may explain why MCA is not the same as the one in the article.
    """
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/awk/extract_autosomes/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
    output:
        "out/chromdet/redo_original_paper_table_s1_autosomes/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1_autosomes"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v

        """


rule chromdet_redo_original_paper_table_s1_autosomes_forced_7_clusters:
    """
    Created:
        2018-04-30 11:32:10
    Aim:
        segmentation files were containing all chromosomes whereas material and methods explicitly use only autosomes. This may explain why MCA is not the same as the one in the article.
    Note:
        Does not work:
                Number of axes selected: 1
                        Percentage of initial variance considered informative: 41.6454 %
                                        Percentage of variance explained by the axis 1  (selected): 41.6454 %
                                                Performing k-means clustering for a fixed number of groups equal to: 7
                                                        No stable cluster results have been found
                                                        RUNNING - Step 2: Extracting S3det results
                                                        No stable clustering was detected

    """
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/awk/extract_autosomes/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
    output:
        "out/chromdet/redo_original_paper_table_s1_autosomes_forced_7_clusters/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1_autosomes_forced_7_clusters"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v -f "-v -k 7"

        """



rule chromdet_redo_original_paper_table_s1_with_disease_autosomes:
    """
    sm out/chromdet/redo_original_paper_table_s1_autosomes/samples_collapsed_filtered_chromatin_space.tsv out/chromdet/redo_original_paper_table_s1_with_disease_autosomes/samples_collapsed_filtered_chromatin_space.tsv
    """
    input:
        dir="out/unzip/d/wget/https/github.com/david-juan/ChromDet/archive/master",
        #segments=expand("out/ln/chromdet_redo_original_paper/{sample}.bed", sample=BASENAME_SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1)
        segments=expand("out/awk/extract_autosomes/tar/xvzf_chromdet_redo_original_paper/{sample}.bed", sample=SEGMENT_SAMPLES_CHROMDET_REDO_ORIGINAL_PAPER_TABLE_S1_WITH_DISEASE)
    output:
        "out/chromdet/redo_original_paper_table_s1_with_disease_autosomes/samples_collapsed_filtered_chromatin_space.tsv"
    params:
        outdir="out/chromdet/redo_original_paper_table_s1_with_disease_autosomes"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        for SAMPLE in {input.segments}
        do
            SAMPLE_BASENAME=`basename $SAMPLE`
            ln -f $SAMPLE {params.outdir}/$SAMPLE_BASENAME
        done

        set +u; source opt/miniconda/bin/activate chromdet; set -u
        export CPLUS_INCLUDE_PATH={WDIR}/opt/miniconda/envs/chromdet/include
        export LD_RUN_PATH={WDIR}/opt/miniconda/envs/chromdet/lib
        export PATH={WDIR}/opt/miniconda/envs/chromdet/bin:/usr/bin

        cd {input.dir}/ChromDet-master/scripts
        ./run_S3det_analysis.pl -d {WDIR}/{params.outdir} -c ../test/States_collapse.txt -s {WDIR}/opt/S3Det_modified/ -v

        """

