"""
Importante note:
    gemBS drastically changed between v2 and v3.
    Below are rules working with gemBS v2.
    gemBS v2 is not available on bioconda so rules below are not reproducible.
    I should spend time to rewrite them.
    http://statgen.cnag.cat/GEMBS/v2/UserGuide/_build/html/pipelineSteps.html
    http://statgen.cnag.cat/GEMBS/UserGuide/_build/html/pipelineIndex.html

Deprecated note:
    Remember to export lib before using bs_call:
    source /gpfs/projects/spicuglia/mw/opt/miniconda/bin/activate gembs
    export LD_LIBRARY_PATH=/gpfs/projects/spicuglia/mw/opt/miniconda/envs/gembs/lib

Blueprint protocol for bisulfite-Seq:
    http://dcc.blueprint-epigenome.eu/#/md/bs_seq_grch38

GRCh38 reference for BS:
    ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
"""

rule gemBS_bscall_r:
    """
    Created:
        2018-04-03 14:06:24
    Important note: 
        The wrapper using gemBS bscall to bs_call does not seems to work correctly.
    Test:
        out/gemBS/bscall_r-GRCh38-BP-BS-protocol_p/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028_chr19.bcf
        """
    input:
        bam = "out/{filler}.bam",
        bai = "out/{filler}.bam.bai",
        fna = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
    output:
        bcf = "out/gemBS/bscall{_p}_r-{fa_genome_id}/{filler}_{chr}.bcf"
    params:
        e=params_gemBS_bscall_e # I should write correspondence between reference "-r" and specie "-e" so user does not have to write both as parameters. 
    wildcard_constraints:
        _p="_p|"
    conda:
        "../envs/gembs.yaml"
    threads:
        1 # Some steps can be parallelized but it seems that the bs calling is not so I would rather keep this pipeline on one thread.
    shell:
        """
        SAMPLE=`basename {output.bcf} | sed 's/_{wildcards.chr}.bcf//'`
        echo "SAMPLE: $SAMPLE"
        OUTDIR=`dirname {output.bcf}`
        echo "OUTDIR: $OUTDIR"

        gemBS\
            bscall\
            --input-bam {input.bam}\
            -r {input.fna}\
            --paired-end\
            -e {params.e}\
            -s $SAMPLE\
            -c {wildcards.chr}\
            --output-dir $OUTDIR
            
        """
      
rule gemBS_bscall_concatenate:
    """
    Created:
        2018-04-03 14:06:24
    Important note: 
        The wrapper using gemBS bscall to bs_call does not seems to work correctly.
    Test:
        out/gemBS/bscall-concatenate_p_r-GRCh38-BP-BS-protocol/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028.bcf
        """
    input:
        #bcf_list should look reference genome and expand the appropriate chromosomes.
        bcf_list=input_gemBS_bscall_concatenate_bcf_list,
        gemBS="opt/miniconda/envs/gembs/bin/gemBS",
        samtools="opt/miniconda/envs/gembs/bin/samtools"
    output:
        bcf="out/gemBS/bscall-concatenate{_p}_r-{fa_genome_id}/{filler}.raw.bcf"
    wildcard_constraints:
        _p="_p|"
    params:
        ld_path="opt/miniconda/envs/gembs/lib",
    threads:
        1 # Some steps can be parallelized but it seems that the bs calling is not so I would rather keep this pipeline on one thread.
    shell:
        """
        set +u; source opt/miniconda/bin/activate gembs; set -u
        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
        
        SAMPLE=`basename {output.bcf}  | sed 's/.raw.bcf$//'`
        echo "SAMPLE: $SAMPLE"
        OUTDIR=`dirname {output.bcf}`
        echo "OUTDIR: $OUTDIR"

        gemBS\
            bscall-concatenate\
            -s $SAMPLE\
            -l {input.bcf_list}\
            -o $OUTDIR
            
        """
 
rule gemBS_methylation_filtering:
    """
    Created:
        2018-04-05 00:39:21
    Test:
        out/gemBS/methylation-filtering/gemBS/bscall-concatenate_p_r-GRCh38-BP-BS-protocol/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028_cpg.txt.gz
        """
    input:
        bcf="out/{filler}.raw.bcf",
        gemBS="opt/miniconda/envs/gembs/bin/gemBS",
        samtools="opt/miniconda/envs/gembs/bin/samtools"
    output:
        bcf="out/gemBS/methylation-filtering/{filler}_cpg.txt.gz"
    wildcard_constraints:
        _p="_p|"
    params:
        ld_path="opt/miniconda/envs/gembs/lib",
    threads:
        1 # Some steps can be parallelized but it seems that the bs calling is not so I would rather keep this pipeline on one thread.
    shell:
        """
        set +u; source opt/miniconda/bin/activate gembs; set -u
        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
        # A LC with '.' for decimal separator should be set here in order to produce 
        export LC_ALL=C
        
        #SAMPLE=`basename {output.bcf} | sed 's/.raw.bcf$//'`
        #echo "SAMPLE: $SAMPLE"
        OUTDIR=`dirname {output.bcf}`
        echo "OUTDIR: $OUTDIR"

        gemBS\
            methylation-filtering\
            -b {input.bcf}\
            -o $OUTDIR
        """


rule gemBS_cpg_bigwig:
    """
    Created:
        2018-04-05 12:50:51
    Test:
        out/gemBS/cpg-bigwig_l-hg38/gemBS/methylation-filtering/gemBS/bscall-concatenate_p_r-GRCh38-BP-BS-protocol/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028_bs_call.bw
        """
    input:
        txt="out/{filler}_cpg.txt.gz",
        gemBS="opt/miniconda/envs/gembs/bin/gemBS",
        samtools="opt/miniconda/envs/gembs/bin/samtools",
        chrominfo = lambda wildcards: config['ids'][wildcards.chrominfo_id]
    output:
        bw_call = "out/gemBS/cpg-bigwig_l-{chrominfo_id}/{filler}.bs_call.bw",
        bw_cov  = "out/gemBS/cpg-bigwig_l-{chrominfo_id}/{filler}.bs_cov.bw"
    params:
        ld_path="opt/miniconda/envs/gembs/lib",
    threads:
        1
    shell:
        """
        set +u; source opt/miniconda/bin/activate gembs; set -u
        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
        export LC_ALL=C
 
        SAMPLE=`basename {output.bw_call} | sed 's/.bs_call.bw$//'`
        #echo "SAMPLE: $SAMPLE"
        OUTDIR=`dirname {output.bw_call}`
        echo "OUTDIR: $OUTDIR"

        gemBS\
            cpg-bigwig\
            -c {input.txt}\
            -l {input.chrominfo}\
            -n $SAMPLE\
            -o $OUTDIR
        """

rule gemBS_test:
    """
    Created:
        2018-04-03 16:50:25
    """
    input:
        tar="opt/miniconda/envs/gembs/gemBS/test/example.tar.gz",
        gemBS="opt/miniconda/envs/gembs/bin/gemBS",
        samtools="opt/miniconda/envs/gembs/bin/samtools"
    output:
        done="out/gemBS/test/done"
    params:
        ld_path="opt/miniconda/envs/gembs/lib",
        outdir="out/gemBS/test"
    shell:
        """
        #Activating gemBS environment
        set +u; source opt/miniconda/bin/activate gembs; set -u
        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
        # A LC with '.' for decimal separator should be set here in order to produce 
        export LC_ALL=C
                
        #Extracting example data
        tar xvzf {input.tar} --directory {params.outdir} --strip 1
        cd {params.outdir}
        touch done

        #Build Index for the reference
        cd reference
        gemBS index -i yeast.fa
        cd ../
    
        #Create JSON metadata file from CSV file
        gemBS prepare-config -t metadata.csv -j metadata.json
    
        #Get list of command mappings to be performed
        gemBS mapping-commands -I reference/yeast.BS.gem -j metadata.json -i ./fastq/ -o ./data/mappings/ -d ./tmp/ -t 8 -p
    
        #Run mapping commands
        gemBS mapping -I reference/yeast.BS.gem -f flowcell_1_1 -j metadata.json -i ./fastq/ -o ./data/mappings/ -d ./tmp/ -t 8 -p
    
        #Merge mappings
        gemBS merging-sample -i ./data/mappings/ -j metadata.json -s test -t 8 -o ./data/sample_mappings/ -d ./tmp/
    
        #Perform Methylation and SNP Calls per Sample and Chromosome
        gemBS bscall -r reference/yeast.fa -e Yeast -s test -c chrIII -i ./data/sample_mappings/test.bam -o ./data/chr_snp_calls/
    
        #Merge Chromosome Calls
        gemBS bscall-concatenate -s test -l ./data/chr_snp_calls/test_chrIII.bcf -o ./data/merged_calls/
    
        #Methylation Filtering
        gemBS methylation-filtering -b data/merged_calls/test.raw.bcf -o ./data/filtered_meth_calls/
    
        #Bisulfite Mapping Report
        gemBS bsMap-report -j metadata.json -i ./data/mappings/ -n TEST_GEMBS -o ./data/report/mappings/
    
        #Variants Report
        gemBS variants-report -l chrIII -j metadata.json -i ./data/chr_snp_calls/ -n TEST_GEMBS -o ./data/report/variants/
    
        #Genome Browser Track
        gemBS cpg-bigwig -c ./data/filtered_meth_calls/test_cpg.txt.gz -l ./reference/chr.len -n TEST_GEMBS -o ./data/Tracks/
    
        """

rule gembs_bs_call_blueprint_like:
    """
    Created:
        2018-04-04 16:21:55
    Important note: 
        The wrapper using gemBS bscall to bs_call does not seems to work correctly.
        Note the difference between bscall (gemBS wrapper) and bs_call (standalone tool).
    Test:
        out/gemBS/bs_call_blueprint_like/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028.bcf
        """
    input:
        bam="out/{filler}.bam",
        bai="out/{filler}.bam.bai",
        fna="out/wget/ftp/ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        bs_call="opt/miniconda/envs/gembs/bin/bs_call",
        samtools="opt/miniconda/envs/gembs/bin/samtools"
    output:
        bcf="out/gemBS/bs_call_blueprint_like/{filler}.bcf"
    params:
        ld_path="opt/miniconda/envs/gembs/lib"
    shell:
        """
        set +u; source opt/miniconda/bin/activate gembs; set -u
        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
        
        samtools\
            view\
            -h\
            {input.bam} |\
            bs_call\
                -r {input.fna}\
                --paired-end |\
                bcftools\
                    convert\
                    -o {output.bcf}\
                    -O b
        """

#rule gembs_bs_call_blueprint_like:
#    """
#    Created:
#        2018-04-04 16:21:55
#    Important note: 
#        The wrapper using gemBS bscall to bs_call does not seems to work correctly.
#        Note the difference between bscall (gemBS wrapper) and bs_call (standalone tool).
#    Test:
#        out/gemBS/bs_call_blueprint_like/inp/bam/GRCh38/Blueprint/TH110_CD34_Bisulfite_NECH0002_1.BS.gem_cnag_bs.GRCh38.20151028.vcf
#        """
#    input:
#        bam="out/{filler}.bam",
#        bai="out/{filler}.bam.bai",
#        fna="out/wget/ftp/ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
#        filter_vcf="opt/miniconda/envs/gembs/bin/filter_vcf",
#        samtools="opt/miniconda/envs/gembs/bin/samtools"
#    output:
#        vcf="out/gemBS/bs_call_blueprint_like/{filler}.vcf"
#    params:
#        ld_path="opt/miniconda/envs/gembs/lib"
#    shell:
#        """
#        set +u; source opt/miniconda/bin/activate gembs; set -u
#        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
#        
#        filter_vcf {input.vcf} 
#        """
#
rule gembs_methylation_filtering:
    """
    Created:
        2018-04-04 17:23:41
    """
    input:
        bcf="out/{filler}.bcf",
        gemBS="opt/miniconda/envs/gembs/bin/gemBS"
    output:
        txt="out/gemBS/methylation_filtering/{filler}.txt.gz"
    params:
        ld_path="opt/miniconda/envs/gembs/lib",
    shell:
        """
        set +u; source opt/miniconda/bin/activate gembs; set -u
        export LD_LIBRARY_PATH={WDIR}/{params.ld_path}
    
        OUTDIR=`dirname {output.txt}`
    
        gemBS methylation-filtering -b {input.bcf} -o $OUTDIR 
        """

