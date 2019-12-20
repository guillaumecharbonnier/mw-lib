rule chipseqspike_test_example:
    """
    Created:
        2018-02-28 14:16:45
    Aim:
        Testing ChIPSeqSpike package.       
    """
    input:
        Rscript="opt/miniconda/envs/r-devel-xml/bin/Rscript",
        script="out/wget/https/bioconductor.org/packages/3.7/bioc/vignettes/ChIPSeqSpike/inst/doc/ChIPSeqSpike.R",
    output:
        "out/chipseqspike/test_example/Rplots.pdf"
    shell:
        """
        cd out/chipseqspike/test_example

        {WDIR}/{input.Rscript} {WDIR}/{input.script}
        """

rule chipseqspike_test_example_docopt:
    """
    Created:
        2018-03-01 13:56:47
    Aim:
        Testing ChIPSeqSpike package, using docopt script.
    """
    input:
        Rscript="opt/miniconda/envs/r-devel-xml/bin/Rscript",
        script="src/r/script/ChIPSeqSpike.R",
        csv="opt/miniconda/envs/r-devel-xml/lib64/R/library/ChIPSeqSpike/extdata/info.csv",
        gff="opt/miniconda/envs/r-devel-xml/lib64/R/library/ChIPSeqSpike/extdata/test_coord.gff"
    params:
        bam="opt/miniconda/envs/r-devel-xml/lib64/R/library/ChIPSeqSpike/extdata/bam_files",
        bw ="opt/miniconda/envs/r-devel-xml/lib64/R/library/ChIPSeqSpike/extdata/bigwig_files",
        out="out/chipseqspike/test_example_docopt",
        assembly="hg19"
    output:
        done="out/chipseqspike/test_example_docopt/Rplots.done"
    shell:
        """
        {input.Rscript}\
            {input.script}\
            -i {input.csv}\
            -b {params.bam}\
            -w {params.bw}\
            -g {input.gff}\
            -a {params.assembly}\
            -o {params.out}

        touch {output.done}
        """

def input_dep_chipseqspike(wildcards):
    """
    Created:
        2018-08-22 16:42:26
    Aim:
    """
    csv=wildcards['csv']

    df = pandas.read_csv("src/chipseqspike/"+csv+".csv", comment='#', delimiter=",")
    #print(df['endogenousBam','exogenousBam','inputBam','bigWigEndogenous','bigWigInput'])
    #paths=df['endogenousBam','exogenousBam','inputBam','bigWigEndogenous','bigWigInput'].tolist()
    paths=df['endogenousBam'].tolist() + df['exogenousBam'].tolist() + df['inputBam'].tolist() + df['bigWigEndogenous'].tolist() + df['bigWigInput'].tolist()

    print(paths)
    return(paths)


rule chipseqspike_gff_csv:
    """
    Created:
        2018-03-01 13:56:47
    Aim:
        Testing ChIPSeqSpike package, using docopt script.
    Note:
        Sam outputs were converted to Bam with Samtools v1.0.6 (Li et al. 2009) and sorted with Picard tools v1.88. Data were
        further processed with Pasha v0.99.21 (Fenouil et al. 2009). Fixed steps wiggle files were
        converted to bigwigs with wigToBigWig.
    Test:
        out/chipseqspike/gff-mm10/test_GRCm38/done
        out/chipseqspike/gff-hg38-Ensembl-93-gene-chipseqspike-example-like/fastkd1_pe_vsl/info.csv
    """
    input:
        Rscript = "opt/miniconda/envs/r-devel-xml/bin/Rscript",
        script = "src/r/script/ChIPSeqSpike.R",
        #bam="",
        #bw="",
        dep = input_dep_chipseqspike,
        csv = "src/chipseqspike/{csv}.csv",
        gff = lambda wildcards: eval(config['ids'][wildcards.gff_id])
    params:
        bam="out/chipseqspike/gff-{gff_id}/{csv}",
        bw ="out/chipseqspike/gff-{gff_id}/{csv}",
        out="out/chipseqspike/gff-{gff_id}/{csv}",
        csv="out/chipseqspike/gff-{gff_id}/{csv}/info.csv"
    output:
        csv="out/chipseqspike/gff-{gff_id}/{csv}/info.csv",
        pdf=expand("out/chipseqspike/gff-{{gff_id}}/{{csv}}/{pdf}.pdf", pdf=[
            #"plotCor_log",
            #"plotCor_number",
            #"plotCor",
            "boxplotSpike",
            "boxplotSpike_violin",
            "plotHeatmap",
            "plotProfile_notScaled",
            "plotProfile",
            "plotTransform"])
    threads:
        MAX_THREADS #using only one but huge memory consumption.
    shell:
        """
        # Extracting genome name from gff_id
        # This require that the first element of gff_id is the genome name as wanted by R package.
        # e.g. hg38-Ensembl-93-gene-chipseqspike-example-like
        A=`echo "{wildcards.gff_id}" | cut -f1 -d'-'`
        
        rm -rf {params.bam} {params.bw}
        mkdir -p {params.bam} {params.bw}
        
        ln {input.script} {params.out}/ChIPSeqSpike.R
        ln {input.gff} {params.out}/anno.gff3

        awk 'NR > 1' {input.csv} > {output.csv}

        for FILE in `cut -f2 -d ',' {output.csv}`;
        do
            echo "$FILE {params.bam}/endo_$(basename $FILE)"
            ln $FILE {params.bam}/endo_$(basename $FILE)
        done

        echo "ENDO done"

        for FILE in `cut -f3 -d ',' {output.csv}`;
        do
            ln $FILE {params.bam}/exo_$(basename $FILE)
        done

        for FILE in `cut -f4 -d ',' {output.csv}`;
        do
            ln -f $FILE {params.bam}/inp_$(basename $FILE)
        done

        for FILE in `cut -f5 -d ',' {output.csv}`;
        do
            ln $FILE {params.bw}/endo_$(basename $FILE)
        done

        for FILE in `cut -f6 -d ',' {output.csv}`;
        do
            ln -f $FILE {params.bw}/inp_$(basename $FILE)
        done


        awk 'BEGIN {{FS=OFS=","}}; {{
            gsub(/out\/.*\//,"endo_",$2);
            gsub(/out\/.*\//,"exo_",$3);
            gsub(/out\/.*\//,"inp_",$4);
            gsub(/out\/.*\//,"endo_",$5);
            gsub(/out\/.*\//,"inp_",$6);

            print}}' {input.csv} > {output.csv}

        echo "output.csv:"
        cat {output.csv}
        #ln {output.csv} {params.csv}

        set +u; source opt/miniconda/bin/activate r-devel-xml; set -u
        cd {params.out}

        echo "Rscript ChIPSeqSpike.R -i info.csv -b '.' -w '.' -g anno.gff3 -a $A -o '.' -p" > script.sh
        bash script.sh
        #    {WDIR}/{input.script}\
        #    -i {WDIR}/{output.csv}\
        #    -b '.'\
        #    -w '.'\
        #    -g {WDIR}/{input.gff}\
        #    -a {wildcards.gff_id}\
        #    -o '.'


        #{WDIR}/{input.Rscript}\
        #    {WDIR}/{input.script}\
        #    -i {WDIR}/{output.csv}\
        #    -b '.'\
        #    -w '.'\
        #    -g {WDIR}/{input.gff}\
        #    -a {wildcards.gff_id}\
        #    -o '.'


        #{input.Rscript}\
        #    {input.script}\
        #    -i {output.csv}\
        #    -b {params.bam}\
        #    -w {params.bw}\
        #    -g {input.gff}\
        #    -a {wildcards.gff_id}\
        #    -o {params.out}

        """





