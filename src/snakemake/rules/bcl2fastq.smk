
#out/bcl2fastq/Run_310_NS500-217_25_02_2020_SK/S003685_1_input_WT/1_input_WT_S1_L001_R1_001.fastq.gz
#
#bcl2fastq_dir = "out/bcl2fastq/"
#output_fastq_prefix = ["S003685_1_input_WT/1_input_WT", "S003686_2_input_G3/2_input_G3", "S003687_3_input_G5/3_input_G5", "S003688_5_ChIP_H4K5Ac_WT/5_ChIP_H4K5Ac_WT", "S003689_5_ChIP_H4K5Ac_G3/5_ChIP_H4K5Ac_G3", "S003690_6_ChIP_H4K5Ac_G5/6_ChIP_H4K5Ac_G5", "S003691_7_ChIP_H4K5Cro_WT/7_ChIP_H4K5Cro_WT", "S003692_8_ChIP_H4K5Cro_G3/8_ChIP_H4K5Cro_G3", "S003693_9_ChIP_H4K5Cro_G5/9_ChIP_H4K5Cro_G5"]
#
#i = 1
#
#int_list = []
#for prefix in output_fastq_prefix:
#    for j in range(1,5):
#        int_list.append(bcl2fastq_dir + "Run_310_NS500-217_25_02_2020_SK/" + prefix + "_S" + str(i) + "_L00" + str(j) + "_R1_001.fastq.gz")
#        int_list.append(bcl2fastq_dir + "Run_310_NS500-217_25_02_2020_SK/" + prefix + "_S" + str(i) + "_L00" + str(j) + "_R2_001.fastq.gz")
#        i = i + 1
#
#def bcl2fastq_input(wildcards):
#    filler = wildcards['filler']
#    accession = 
#
#    paths="out/bcl2fastq/"+filler+"/Reports/html/tree.html"
#    return(paths)
#

rule bcl2fastq:
    """
    Created:
        2020-05-18 00:37:53
    Aim:
        Produce fastq from bcl
        Note the fastq are not explicitely defined here so you should first have a
        Snakemake workflow requiring the tree.html file for all your runs, which will also produce all fastq files.
        Then use a second workflow that knows where the produced fastq are.
        Please open a github issue if you know how to properly explicitely define fastq in this case.
        Maybe checkpoints ?

    Test:
        out/bcl2fastq/_--barcode-mismatches_1/gpfs/projects/spicuglia/mw/inp/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/Reports/html/tree.html out/bcl2fastq/_--barcode-mismatches_1_--no-lane-splitting/gpfs/projects/spicuglia/mw/inp/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/Reports/html/tree.html
    """
    input:
        xml="/{filler}/RunInfo.xml", 
        csv="out/{tool}{extra}/{filler}/SampleSheet.csv",
        #csv="/{filler}/SampleSheet.csv"
        #"inp/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/RunInfo.xml"
        #bcl2fastq_input
    output:
        #[x.strip() for x in open("/"+wildcards["filler"]+"/output_files.txt","r")]
        #directory("out/bcl2fastq/{filler}")
        #"out/bcl2fastq/{accession}/{experiment}/
        #clustering/{sample}")
        xml="out/{tool}{extra}/{filler}/RunInfo.xml",
        html="out/{tool}{extra}/{filler}/Reports/html/tree.html"
        #int_list
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bcl2fastq/"
    threads:
        8
    conda:
        "../envs/bcl2fastq.yml"
    log:
        #bcl2fastq_dir + "bcl2fastq.log"
        bcl2fastq_log = "out/{tool}{extra}/{filler}/log",
        index_in_undetermined = "out/{tool}{extra}/{filler}/indexes.txt"
    shell:
        """
        cp {input.xml} {output.xml}
        INDIR=`dirname {input.xml}`
        OUTDIR=`dirname {output.xml}`
        (bcl2fastq --input-dir $INDIR/Data/Intensities/BaseCalls --runfolder-dir $OUTDIR --output-dir $OUTDIR {params.extra}) &> {log.bcl2fastq_log} 
        (zcat $OUTDIR/Undetermined_S0_R1_001.fastq.gz | grep '^@' | cut -d : -f 10 | sort | uniq -c | sort -nr > $OUTDIR/indexes.txt ) &> {log.index_in_undetermined}
        """

#find . -type f -name '*.fastq.gz' -mindepth 2 -exec ln -sf -- {} . \;
#--sample-sheet /gpfs/tgml/reads/analysis_results/202002007/samplesheet.csv \

#rule expose_fastq_from_bcl2fastq:
#    """
#    out/bcl2fastq/gpfs/projects/spicuglia/mw/inp/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/Run_310_NS500-217_25_02_2020_SK/S003693_9_ChIP_H4K5Cro_G5/9_ChIP_H4K5Cro_G5_S9_L001_R1_001.fastq.gz
#    """
#    input:
#        "out/bcl2fastq/{filler}/Reports/html/tree.html"
#    output:
#        "out/bcl2fastq/{filler}/{path}.fastq.gz"
#    shell:
#        "touch {output}"
#
# out/bcl2fastq/Reports/html/tree.html

