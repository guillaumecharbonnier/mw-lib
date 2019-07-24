#rule ln_extra:
#    """
#    Created:
#        2018-11-07 01:54:53
#        tectonic/ln/_-srf/src/tex/
#    """
#    input:
#        "{filler}"
#    output:
#        "out/{tool}{extra}/{filler}"
#    wildcard_constraints:
#        tool="ln/"
#    params:
#        extra = params_extra
#    conda:
#        "../envs/coreutils.yaml"
#    shell:
#        "ln {extra} {input} {output}"
#
localrules: ln_alias
rule ln_alias:
    """
    Created:
        2017-10-03 09:09:59
    Modified:
        2018-08-20 18:00:29 - Adding quotes to allow selection of files with spaces.
        2019-01-28 10:18:00 - Changed to symbolic link instead.
    Test:
        out/ln/alias/experiments/encode_broad_h3k4me3/Normal_Osteo.bed
        out/ln/alias/microarray-ko-nut-log-fc-neg.bed
    """
    input:
        lambda wildcards: config['ids'][wildcards.id]
        #input_ln_alias
    output:
        alias="out/ln/alias/{id}"
    #wildcard_constraints:
    #    alias_id=".*"
    conda:
        "../envs/coreutils.yaml"
    shell:
        # quotes around {input} is mandatory to be able
        # to have files with special characters in input.
        "ln -srf '{input}' {output.alias}"

localrules: ln_dynamic
rule ln_dynamic:
    """
    Created:
        2018-11-12 17:32:11
    Modified:
        2018-11-22 14:24:31 - dirdynsep (directory-dynamic-separator) is added in order to allow Snakemake to differentiate between the static directory and possible subdirectories in dynamic part. This is useful for example for iGenomes.
    Aim:
        Work like "ln_alias" rule but explicitely look for "done" file to know input files are available, as dynamic file are not declared in upstream rule.
        This rule allows to explicitely declare dynamic files produced by a rule which are wanted for downstream rules.
        This is useful with ChromHMM when some plots and bed filenames depend on config.
    Test:
        out/ln/dynamic/ChromHMM/LearnModel_test19_numstates-10_assembly-mm10/meta_10_segments.bed
        {dynamic}", dynamic = ["meta_10_segments.bed","meta_10_dense.bed", "meta_10_expanded.bed", "meta_10_overlap.txt"])
        expand("out/ln/dynamic/ChromHMM/LearnModel_test19_numstates-10_assembly-mm10/{dynamic}", dynamic = ["meta_10_segments.bed","meta_10_dense.bed", "meta_10_expanded.bed", "meta_10_overlap.txt"])

        out/ln/dynamic/tar/xvzf/wget/ftp/igenome:G3nom3s4u@ussd-ftp.illumina.com/Drosophila_melanogaster/UCSC/dm6/Drosophila_melanogaster_UCSC_dm6/dirdynsep/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome.1.bt2
    """
    input:
        directory = "out/{dir}/done"
    params:
        dynamic = "out/{dir}/{dynamic}"
    output:
        alias = "out/ln/dynamic/{dir}/dirdynsep/{dynamic}"
    conda:
        "../envs/coreutils.yaml"
    shell:
        #"ln -L {params.dynamic} {output.alias}"
        "ln -srf {params.dynamic} {output.alias}"

localrules: ln_srf_project_dir
rule ln_srf_project_dir:
    """
    Created:
        2018-10-14 03:06:07
    Aim:
       
    """
    input:
        "{filler}"
    output:
        "out/ln/projdir/{filler}"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -srf '{input}' '{output}'"

localrules: ln_srf_parent_dir
rule ln_srf_parent_dir:
    """
    Created:
        2019-01-22 15:05:48
    Aim:
        Give access to paths relative to the parent directory.
        |_mw
        |_mw-legacy
    Test:
        out/ln/updir/mw-legacy/inp/testln.tmp
    """
    input:
        "../{filler}"
    output:
        "out/ln/updir/{filler}"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -srf '{input}' '{output}'"

localrules: ln_srf_abspath
rule ln_srf_abspath:
    """
    Created:
        2019-01-22 15:05:48
    Aim:
    Test:
        out/ln/abspath/gpfs/projects/spicuglia/mw-legacy/inp/testln.tmp
    """
    input:
        "/{filler}"
    output:
        "out/ln/abspath/{filler}"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -srf '{input}' '{output}'"

localrules: ln_broadPeak_to_bed
rule ln_broadPeak_to_bed:
    """
    Created:
        2018-11-23 15:00:41
    Aim:
        For most usages, broadPeak are like bed so here is a simple alias to be able to use all bed related rules.
    Test:
        out/ln/broadPeak_to_bed/gunzip/to-stdout/wget/ftp/hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneOsteoH3k04me3Pk.bed
    """
    input:
        "out/{filler}.broadPeak"
    output:
        "out/ln/broadPeak_to_bed/{filler}.bed"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -srf {input} {output}"

localrules: ln_renamed_blueprint_fastq
rule ln_renamed_blueprint_fastq:
    """
    Created:
        2018-04-17 21:32:28
    Aim:
        Rename blueprint fastq according to Salva's need.
    Test:
        expand("out/ln/renamed_blueprint_fastq/{alias_id}", alias_id=BLUEPRINT_ALIAS_SALVA)
    """
    input:
        input_ln_renamed_blueprint_fastq
    output:
        alias="out/ln/renamed_blueprint_fastq/{alias_id}"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -srf {input} {output.alias}"

localrules: ln_pe_remove_mate_prefix
rule ln_pe_remove_mate_prefix:
    """
    Created:
        2017-12-08 12:18:35
    Aim:
        Sometimes paired-end mates are identified with suffix 'R1' or just '1'. This is very annoying to handle boths until alignment so this rule create a link.
    """
    input:
        fastq="out/{filler}_R{mate}.{ext}"
    output:
        fastq="out/ln/pe_remove_mate_prefix/{filler}_{mate}.{ext}"
    wildcard_constraints:
        mate="1|2",
        ext="fastq|fastq.gz"
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -srf {input.fastq} {output.fastq}"

