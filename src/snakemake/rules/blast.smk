rule blastn:
    """
    Created:
        2017-04-13 11:22:37
    Modified:
        2017-04-26 14:36:46 - Validated feature: strand 'plus' and 'minus' instead of default 'both' in order to see strand specific transcription in repeats.
    Note:
        Blast best hit only:
            http://seqanswers.com/forums/showthread.php?t=23166
            blastn -query transcripts.fa -out transcripts.blast.txt -task megablast -db refseq_rna -num_threads 12 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 7 -perc_identity 50 -max_target_seqs 1 &
            Adding the "-max_target_seqs" flag and setting it to "1" yields what appears to be the best hit in terms of e-value and bit score. I haven't done extensive comparisons, but it appears where multiple matches yield the same e-value (e.g., 0), the match with the highest score is retained.
    Doc:
        blastn [-h] [-help] [-import_search_strategy filename]
            [-export_search_strategy filename] [-task task_name] [-db database_name]
            [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
            [-negative_gilist filename] [-entrez_query entrez_query]
            [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
            [-subject subject_input_file] [-subject_loc range] [-query input_file]
            [-out output_file] [-evalue evalue] [-word_size int_value]
            [-gapopen open_penalty] [-gapextend extend_penalty]
            [-perc_identity float_value] [-qcov_hsp_perc float_value]
            [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
            [-xdrop_gap_final float_value] [-searchsp int_value]
            [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
            [-min_raw_gapped_score int_value] [-template_type type]
            [-template_length int_value] [-dust DUST_options]
            [-filtering_db filtering_database]
            [-window_masker_taxid window_masker_taxid]
            [-window_masker_db window_masker_db] [-soft_masking soft_masking]
            [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
            [-best_hit_score_edge float_value] [-window_size int_value]
            [-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
            [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
            [-outfmt format] [-show_gis] [-num_descriptions int_value]
            [-num_alignments int_value] [-line_length line_length] [-html]
            [-max_target_seqs num_sequences] [-num_threads int_value] [-remote]
            [-version]
    Test:
        out/blast/blastn_task-blastn_db-repbase-mus-musculus_evalue-1e-3_word_size-7_best_hit_score_edge-0.05_best_hit_overhang-0.25_perc_identity-50_strand-plus_outfmt-1_max_target_seqs-1/seqtk/seq_A/star/pe_GRCm38_outFilterMultimapNmax-1/sickle/pe_-t_sanger_-q_20/merge_lanes/run176/RNA-S-H2AL2-WT-Rep1/unmapped_pair1.tsv
    """
    input:
        db_files=expand("out/blast/makeblastdb/Mus_musculus_all.{ext}", ext=["nhr","nin","nog","nsd","nsi","nsq"]),
        fasta="out/{filler}.fasta"
    output:
        tsv="out/blast/blastn_task-{task}_db-repbase-mus-musculus_evalue-{evalue}_word_size-{word_size}_best_hit_score_edge-{best_hit_score_edge}_best_hit_overhang-{best_hit_overhang}_perc_identity-{perc_identity}_strand-{strand}_outfmt-{outfmt}_max_target_seqs-{max_target_seqs}/{filler}.tsv"
    params:
        db_prefix="out/blast/makeblastdb/Mus_musculus_all"
    wildcard_constraints:
        task="blastn|blastn-short|dc-megablast|megablast|rmblastn",
        evalue="1e-[0-9]+",
        word_size="[0-9]+",
        best_hit_score_edge="0.[0-9]+",
        best_hit_overhang="0.[0-9]+",
        perc_identity="[0-9]+",
        max_target_seqs="[0-9]+",
        strand="plus|minus|both",
        outfmt="[0-9]+"
    threads:
        16
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn \
            -task {wildcards.task} \
            -query {input.fasta} \
            -db {params.db_prefix} \
            -num_threads {threads} \
            -out {output.tsv} \
            -evalue {wildcards.evalue} \
            -word_size {wildcards.word_size} \
            -best_hit_score_edge {wildcards.best_hit_score_edge} \
            -best_hit_overhang {wildcards.best_hit_overhang} \
            -perc_identity {wildcards.perc_identity} \
            -outfmt {wildcards.outfmt} \
            -strand {wildcards.strand} \
            -max_target_seqs {wildcards.max_target_seqs}
        """

rule makeblastdb:
    """
    Created:
        2017-04-13 11:22:31
    """
    input:
        fasta="inp/repbase/Mus_musculus_all.fasta"
    output:
        db_files=expand("out/blast/makeblastdb/Mus_musculus_all.{ext}", ext=["nhr","nin","nog","nsd","nsi","nsq"]) 
    conda:
        "../envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input.fasta} -parse_seqids -dbtype nucl -out out/blast/makeblastdb/Mus_musculus_all
        """
