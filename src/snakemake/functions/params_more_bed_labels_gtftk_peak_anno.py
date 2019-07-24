def params_more_bed_labels_gtftk_peak_anno(wildcards):
    """
    Created:
        2017-05-09 17:32:02
    """
    bed_anno_list_id=wildcards['bed_anno_list_id']

    if bed_anno_list_id == 'mm10-main-chr-tss-classes-and-brd-peaks':
# Test without condition.
#elif bed_anno_list == 'tss_classes_with_without_brdt_brd4_and_h4k5ac':
        #labels="TSS +Brdt +H4K5ac,TSS +Brd4 +H4K5ac,TSS -Brdt -Brd4 +H4K5ac,Brdt peaks,Brd4 peaks"
        labels="TSS_Brdt_H4K5ac,TSS_Brd4_H4K5ac,TSS_no_Brdt_no_Brd4_H4K5ac,Brdt_peaks,Brd4_peaks"
   
    elif bed_anno_list_id == 'mm10-main-chr-nuc-r-cpg-h2alap1-peaks':
        labels="Nuc. pos. in R,CpG islands,H2A.Lap1 peaks"

    elif (bed_anno_list_id == 'gsat-mm-synrep-mm') or (bed_anno_list_id == 'mm10-main-chr-gsat-mm-synrep-mm') or (bed_anno_list_id == 'GRCm38-gsat-mm-synrep-mm'):
        labels="Major satellites,Minor satellites"

    elif bed_anno_list_id == 'mm10-main-chr-all-2017-06-26':
        #labels="TSS +Brdt +H4K5ac,TSS +Brd4 +H4K5ac,TSS -Brdt -Brd4 +H4K5ac,Brdt peaks,Brd4 peaks,Nuc. pos. in R,CpG islands,H2A.Lap1 peaks,Major satellites,Minor satellites"
        labels="TSS_Brdt_H4K5ac,TSS_Brd4_H4K5ac,TSS_no_Brdt_no_Brd4_H4K5ac,Brdt_peaks,Brd4_peaks,Nuc_pos_in_R,CpG_islands,H2ALap1_peaks,Major_satellites,Minor_satellites"

    elif bed_anno_list_id == 'mm10-main-chr-all-2017-09-28':
        #labels="TSS +Brdt +H4K5ac,TSS +Brd4 +H4K5ac,TSS -Brdt -Brd4 +H4K5ac,Brdt peaks,Brd4 peaks,Nuc. pos. in R,CpG islands,H2A.Lap1 peaks,Major satellites,Minor satellites,upreg in Nut-KO,downreg in Nut-KO,TSS+-2kb upreg,TSS+-2kb downreg"
        labels="TSS_Brdt_H4K5ac,TSS_Brd4_H4K5ac,TSS_no_Brdt_no_Brd4_H4K5ac,Brdt_peaks,Brd4_peaks,Nuc_pos_in_R,CpG_islands,H2ALap1_peaks,Major_satellites,Minor_satellites,upreg_in_Nut_KO,downreg_in_Nut_KO,TSSpm2kb_upreg,TSSpm2kb_downreg"

    elif bed_anno_list_id == 'mm10-main-chr-all-2017-10-04':
        #labels="TSS +Brdt +H4K5ac,TSS +Brd4 +H4K5ac,TSS -Brdt -Brd4 +H4K5ac,Brdt peaks,Brd4 peaks,Nuc. pos. in R,CpG islands,H2A.Lap1 peaks,Major satellites,Minor satellites,upreg in Nut-KO,downreg in Nut-KO,TSS+-3kb upreg RNA-Seq,TSS+-3kb downreg RNA-Seq,TSS+-3kb upreg Microarray,TSS+-3kb downreg Microarray"
        labels="TSS_Brdt_H4K5ac,TSS_Brd4_H4K5ac,TSS_no_Brdt_no_Brd4_H4K5ac,Brdt_peaks,Brd4_peaks,Nuc_pos_in_R,CpG_islands,H2ALap1_peaks,Major_satellites,Minor_satellites,upreg_in_Nut_KO,downreg_in_Nut_KO,TSSpm3kb_upreg_RNA_Seq,TSSpm3kb_downreg_RNA_Seq,TSSpm3kb_upreg_Microarray,TSSpm3kb_downreg_Microarray"

    else:
        print('Error: Add labels function params_more_bed_labels_gtftk_peak_anno')

    return(labels)

