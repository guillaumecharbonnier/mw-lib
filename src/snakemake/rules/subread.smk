rule subread_buildindex_extra:
    """
    Created:
        2018-10-24 14:28:34
    Doc:
        http://gensoft.pasteur.fr/docs/subread/1.4.6-p3/SubreadUsersGuide.pdf#page=14
Usage:

 ./subread-buildindex [options] -o <basename> {FASTA file1} [FASTA file2] ...

Required arguments:

    -o <basename>   base name of the index to be created

Optional arguments:

    -F              build a full index for the reference genome. 16bp subreads
                    will be extracted from every position of the reference
                    genome. Size of the index is typically 3 times the size of
                    index built from using the default setting.

    -B              create one block of index. The built index will not be split
                    into multiple pieces. This makes the largest amount of
                    memory be requested when running alignments, but it enables
                    the maximum mapping speed to be achieved. This option
                    overrides -M when it is provided as well.

    -M <int>        size of requested memory(RAM) in megabytes, 8000 by default.

    -f <int>        specify the threshold for removing uninformative subreads
                    (highly repetitive 16mers in the reference). 100 by default.

    -c              build a color-space index.

    -v              output version of the program.
    Note:
        Subread does not accept fasta bigger than 4GB hence the need for a list of fasta for each chr as input...
    Test:
        out/subread/buildindex/fa-genome-GRCh38-r95-chr22.reads
        out/subread/buildindex_-c/fa-genome-GRCh38-r95-chr22.reads
    """
    input:
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id])
        #fai   = lambda wildcards: [path + '.fai' for path in eval(config['ids'][wildcards.bam_list_id])],
        #genome="out/{filler}.fa"
    output:
        expand("out/{{tool}}{{extra}}/{{fa_genome_id}}.{ext}", ext=["reads","log","files","00.c.tab","00.c.array"])
    params:
        o = "out/subread/buildindex{extra}/{fa_genome_id}",
        extra = params_extra
    wildcard_constraints:
        tool="subread/buildindex"
    conda:
        "../envs/subread.yaml"
    threads:
        1
    shell:
        "subread-buildindex {params.extra} -o {params.o} {input.fasta}"

rule subread_align_extra:
    """
    Created:
        2018-10-24 16:07:12
    Doc:
        http://gensoft.pasteur.fr/docs/subread/1.4.6-p3/SubreadUsersGuide.pdf#page=17
    Usage:
./subread-align [options] -i <index_name> -r <input> -t <type> -o <output>
## Mandatory arguments:
  -i <string>       Base name of the index.

  -r <string>       Name of an input read file. If paired-end, this should be
                    the first read file (typically containing "R1"in the file
                    name) and the second should be provided via "-R".
                    Acceptable formats include gzipped FASTQ, FASTQ and FASTA.
                    These formats are identified automatically.

  -t <int>          Type of input sequencing data. Its values include
                      0: RNA-seq data.
                      1: genomic DNA-seq data.

## Optional arguments:
# input reads and output
  
  -o <string>       Name of an output file. By default, the output is in BAM
                    format. Omitting this option makes the output be written to
                    STDOUT.
  
  -R <string>       Name of the second read file in paired-end data (typically
                    containing "R2" the file name).
  
  --SAMinput        Input reads are in SAM format.

  --BAMinput        Input reads are in BAM format.

  --SAMoutput       Save mapping results in SAM format.

# Phred offset
  
  -P <3:6>          Offset value added to the Phred quality score of each read
                    base. '3' for phred+33 and '6' for phred+64. '3' by default.

# thresholds for mapping
  
  -n <int>          Number of selected subreads, 10 by default.

  -m <int>          Consensus threshold for reporting a hit (minimal number of
                    subreads that map in consensus) . If paired-end, this gives
                    the consensus threshold for the anchor read (anchor read
                    receives more votes than the other read in the same pair).
                    3 by default

  -p <int>          Consensus threshold for the non- anchor read in a pair. 1 by
                    default.

  -M <int>          Maximum number of mis-matched bases allowed in each reported
                    alignment. 3 by default. Mis-matched bases found in soft-
                    clipped bases are not counted.

# unique mapping and multi-mapping
                                                    
  --multiMapping    Report multi-mapping reads in addition to uniquely mapped
                    reads. Use "-B" to set the maximum number of equally-best
                    alignments to be reported.

  -B <int>          Maximum number of equally-best alignments to be reported for
                    a multi-mapping read. Equally-best alignments have the same
                    number of mis-matched bases. 1 by default.

# indel detection

  -I <int>          Maximum length (in bp) of indels that can be detected. 5 by
                    default. Indels of up to 200bp long can be detected.

  --complexIndels   Detect multiple short indels that are in close proximity
                    (they can be as close as 1bp apart from each other).

# read trimming
  
  --trim5 <int>     Trim off <int> number of bases from 5' end of each read. 0
                    by default.
  
  --trim3 <int>     Trim off <int> number of bases from 3' end of each read. 0
                    by default.

# distance and orientation of paired end reads

  -d <int>          Minimum fragment/insert length, 50bp by default.

  -D <int>          Maximum fragment/insert length, 600bp by default.

  -S <ff:fr:rf>     Orientation of first and second reads, 'fr' by default (
                    forward/reverse).

# number of CPU threads

  -T <int>          Number of CPU threads used, 1 by default.

# read group

  --rg-id <string>  Add read group ID to the output.

  --rg <string>     Add <tag:value> to the read group (RG) header in the output.

# color space reads

  -b                Convert color-space read bases to base-space read bases in
                    the mapping output. Note that read mapping is performed at
                    color-space.

# dynamic programming

  --DPGapOpen <int> Penalty for gap opening in short indel detection. -1 by
                    default.

  --DPGapExt <int>  Penalty for gap extension in short indel detection. 0 by
                    default.

  --DPMismatch <int> Penalty for mismatches in short indel detection. 0 by
                    default.

  --DPMatch <int>   Score for matched bases in short indel detection. 2 by
                    default.

# detect structural variants

  --sv              Detect structural variants (eg. long indel, inversion,
                    duplication and translocation) and report breakpoints. Refer
                    to Users Guide for breakpoint reporting.

# gene annotation

  -a                Name of an annotation file. GTF/GFF format by default. See
                    -F option for more format information.

  -F                Specify format of the provided annotation file. Acceptable
                    formats include 'GTF' (or compatible GFF format) and
                    'SAF'. 'GTF' by default. For SAF format, please refer to
                    Users Guide.

  -A                Provide a chromosome name alias file to match chr names in
                    annotation with those in the reads. This should be a two-
                    column comma-delimited text file. Its first column should
                    include chr names in the annotation and its second column
                    should include chr names in the index. Chr names are case
                    sensitive. No column header should be included in the
                    file.

  --gtfFeature <string>  Specify feature type in GTF annotation. 'exon'
                    by default. Features used for read counting will be 
                    extracted from annotation using the provided value.

  --gtfAttr <string>     Specify attribute type in GTF annotation. 'gene_id'
                    by default. Meta-features used for read counting will be 
                    extracted from annotation using the provided value.

    Test:
        out/subread/buildindex_-c/GRCh38-r94-chr1.done
        out/subread/buildindex_-c/GRCh38-r94-main-chr.done
    """
    input:
        fasta = lambda wildcards: eval(config['ids'][wildcards.fa_genome_id]),
        index = lambda wildcards: eval(config['ids'][wildcards.subread_index_id]),
        gtf   = lambda wildcards: eval(config['ids'][wildcards.gtf_id])
        #genome="out/{filler}.fa"
    output:
        done="out/{tool}{extra}_-i_{index_id}_-t_{t}/{fa_genome_id}.done"
    params:
        #index_basename = params_subread_index_basename,
        #"out/subread/buildindex{extra}/{fa_genome_id}",
        extra = params_extra
    wildcard_constraints:
        tool="subread/align",
        t="0|1" #0: RNA-seq data, #1:genomic DNA-seq data.
    conda:
        "../envs/subread.yaml"
    threads:
        1
    shell:
        """
        subread-align {params.extra} -i {parmas.index_basename} TODO WHEN CUTADAPT IS DONE.                                                                                                                                                                                              
  -r <string>       Name of an input read file. If paired-end, this should be                                                                                                                 
                    the first read file (typically containing "R1"in the file                                                                                                                 
                    name) and the second should be provided via "-R".                                                                                                                         
                    Acceptable formats include gzipped FASTQ, FASTQ and FASTA.                                                                                                                
                    These formats are identified automatically.                                                                                                                               
                                                                                                                                                                                              
  -t <int>          Type of input sequencing data. Its values include                                                                                                                         
                      0: RNA-seq data                                                                                                                                                         
                      1: genomic DNA-seq data.      -o {params.o} {input.fasta}
        touch {output.done}
        """

