rule sh_generate_scipio_white_list:
    """
    Aim:
        Script from Scipio-bioscience_Cytonaut-Local-UserGuide_v1.3.pdf to generate whitelist for STARsolo
    Test:
        out/sh/generate_scipio_white_list/ln/updir/mw-sst/inp/Raw_data/Scipio/run3_565/MLP_B_S1_R1_L001.white_list.tsv
    """
    input:
        fastq = "out/{filler}.fastq.gz"
    output:
        white_list_count = "out/sh/generate_scipio_white_list/{filler}.white_list_count.tsv",
        white_list = "out/sh/generate_scipio_white_list/{filler}.white_list.tsv"
    shell:
        """
        L_FASTQ="{input.fastq}"
        PATTERN="CCCCCCCCCCCCNNNNNNNNNNNNN"
        output="{output.white_list_count}"
        solo_output="{output.white_list}"
        echo "Script to convert fastq R1 at ${{L_FASTQ}} into white-list-total with optional downsampling."
        echo "This version assumes no downsampling."
        echo "Determining barcode and UMI start and end position in R1s."
        barcode_start=`echo ${{PATTERN}} | grep -b -o C | sed 's/:.*$//' | head -n 1 | sed 's/$/+1/' | bc -l`
        barcode_end=`echo ${{PATTERN}} | grep -b -o C | sed 's/:.*$//' | tail -n 1 | sed 's/$/+1/' | bc -l`
        umi_start=`echo ${{PATTERN}} | grep -b -o N | sed 's/:.*$//' | head -n 1 | sed 's/$/+1/' | bc -l`
        umi_end=`echo ${{PATTERN}} | grep -b -o N | sed 's/:.*$//' | tail -n 1 | sed 's/$/+1/' | bc -l`
        echo "Sample being converted to whitelist:"
        echo ${{L_FASTQ}}
        zcat ${{L_FASTQ}} | \
            awk '(NR%4==2)' | \
            cut -c ${{barcode_start}}-${{barcode_end}},${{umi_start}}-${{umi_end}} --output-delimiter '_' | \
            awk '!seen[$0]++' | \
            cut -d'_' -f1 | \
            awk -F, '{{a[$1]++;}}END{{for (i in a)print i, a[i];}}' | \
            sort -k 2 -n -r > $output
        echo "Generate solo output:"
        echo $solo_output
        cut -d" " -f1 $output > $solo_output
        """