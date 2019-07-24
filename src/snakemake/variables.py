# Acronym for danpos output.
PPR = ["positions","peaks","regions"]

LANES=["L1","L2","L3","L4"]
HUMAN_AUTOSOMES = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
HUMAN_MAIN_CHROMOSOMES = HUMAN_AUTOSOMES + ["X","Y","M"]
HUMAN_MAIN_CHROMOSOMES_NCBI = expand("chr{chr}", chr=HUMAN_MAIN_CHROMOSOMES)
HUMAN_MAIN_CHR_ENSEMBL=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]
MOUSE_MAIN_CHR_ENSEMBL=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y","MT"]

MOUSE_MAIN_CHR=[
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chrX",
    "chrY",
    "chrM"]

TAR_CONTENT_HMCAN_DIFF_TEST_EXAMPLE=[
    "control_C1.sam",
    "chip_C2_rep_3.sam",
    "chip_C1_rep_3.sam",
    "C1_files.txt",
    "C2_files.txt",
    "chip_C1_rep_2.sam",
    "C2_control.txt",
    "chip_C2_rep_2.sam",
    "C1_control.txt",
    "chip_C1_rep_1.sam",
    "chip_C2_rep_1.sam",
    "reference/chr1.fa",
    "data/hg19-blacklist.bed",
    "data/GC_profile_100KbWindow_Mapp76_hg19.cnp",
    "control_C2.sam"]
