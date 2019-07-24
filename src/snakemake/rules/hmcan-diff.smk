rule hmcan_diff_test_example:
    """
    Created:
        2018-02-26 12:08:11
    Aim:
        Test example from: 
        http://www.cbrc.kaust.edu.sa/hmcan/hmcan-diff_desc.php
    """
    input:
        bin="opt/miniconda/envs/hmcan-diff/bin/HMCan-diff",
        files=expand("out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/{filler}", filler=TAR_CONTENT_HMCAN_DIFF_TEST_EXAMPLE)
        #C1_ChIP =       "out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/C1_files.txt",
        #C2_ChIP =       "out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/C2_files.txt",
        #C1_Control =    "out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/C1_control.txt",
        #C2_Control =    "out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/C2_control.txt",
        #GCIndex =       "out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/data/GC_profile_100KbWindow_Mapp76_hg19.cnp",
        #blackListFile = "out/tar/xvzf_hmcan_diff_test_example/hmcan-diff_example/data/hg19-blacklist.bed"
    params:
        genomePath="hmcan-diff_example/reference/",
    shell:
        """
        set +u; source opt/miniconda/bin/activate hmcan-diff; set -u
        
        cd out/tar/xvzf_hmcan_diff_test_example

        HMCan-diff\
            --name hmcan-diff_example\
            --C1_ChIP hmcan-diff_example/C1_files.txt\
            --C2_ChIP hmcan-diff_example/C2_files.txt\
            --C1_Control hmcan-diff_example/C1_control.txt\
            --C2_Control hmcan-diff_example/C2_control.txt\
            --format SAM\
            --genomePath hmcan-diff_example/reference/\
            --GCProfile hmcan-diff_example/data/GC_profile_100KbWindow_Mapp76_hg19.cnp\
            --C1_minLength 145\
            --C1_medLength 150\
            --C1_maxLength 155\
            --C2_minLength 145\
            --C2_medLength 150\
            --C2_maxLength 155\
            --blackListFile hmcan-diff_example/data/hg19-blacklist.bed

        """
