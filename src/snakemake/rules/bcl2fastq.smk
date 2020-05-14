import os

os.system("source activate bcl2fastq")

output_fastq_prefix = ["S003685_1_input_WT/1_input_WT", "S003686_2_input_G3/2_input_G3", "S003687_3_input_G5/3_input_G5", "S003688_5_ChIP_H4K5Ac_WT/5_ChIP_H4K5Ac_WT", "S003689_5_ChIP_H4K5Ac_G3/5_ChIP_H4K5Ac_G3", "S003690_6_ChIP_H4K5Ac_G5/6_ChIP_H4K5Ac_G5", "S003691_7_ChIP_H4K5Cro_WT/7_ChIP_H4K5Cro_WT", "S003692_8_ChIP_H4K5Cro_G3/8_ChIP_H4K5Cro_G3", "S003693_9_ChIP_H4K5Cro_G5/9_ChIP_H4K5Cro_G5"]

i = 1

int_list = []
for prefix in output_fastq_prefix:
	for j in range(1,5):
		int_list.append("/gpfs/tgml/reads/analysis_results/202002007/Run_310_NS500-217_25_02_2020_SK/" + prefix + "_S" + str(i) + "_L00" + str(j) + "_R1_001.fastq.gz")
		int_list.append("/gpfs/tgml/reads/analysis_results/202002007/Run_310_NS500-217_25_02_2020_SK/" + prefix + "_S" + str(i) + "_L00" + str(j) + "_R2_001.fastq.gz")
	i = i + 1

rule bcl2fastq:
	input:               
		"/gpfs/tgml/reads/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/RunInfo.xml"
	output:                           
		int_list
	threads:
		1
	conda:                                                                                        
		"/gpfs/tgml/mw-lib/src/snakemake/envs/bcl2fastq.yml"
	log:                                                                                          
		"/gpfs/tgml/reads/analysis_results/202002007/bcl2fastq.log"
	shell:"""                                                                                    
		INDIR=`dirname {input}`          
		echo $INDIR 
		cp /gpfs/tgml/reads/analysis_results/202002007/samplesheet.csv $INDIR/SampleSheet.csv
		bcl2fastq --input-dir $INDIR/Data/Intensities/BaseCalls --runfolder-dir $INDIR --output-dir /gpfs/tgml/reads/analysis_results/202002007 --barcode-mismatches 1 &> {log}
        """  
#--sample-sheet /gpfs/tgml/reads/analysis_results/202002007/samplesheet.csv \
