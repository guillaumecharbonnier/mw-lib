# rule meme_sea_diff:
#     """
#     Doc:
#         https://meme-suite.org/meme/doc/meme-chip.html?man_type=web
#         http://meme-suite.org/doc/chip.html
#     Warning:
#         If the fasta input is a symlink, the tool will fail to appropriately symlink it again in its output directory...
#     """
#     input:
#         posfa = "out/{filler}/{pos}.fa",
#         negfa   = "out/{filler}/{neg}.fa",
#         meme = lambda wildcards: mwconf['ids'][wildcards.meme_id]
#     output:
#         expand(
#             "out/meme/chip-diff_{{meme_id}}/{{filler}}/{{pos}}_neg_{{neg}}/{file}",
#             file=[
#                 'meme-chip.html',
#                 'summary.tsv'
#             ]
#         )
#     log:
#                "out/meme/chip-diff_{meme_id}/{filler}/{pos}_neg_{neg}/log"
#     benchmark:
#                "out/meme/chip-diff_{meme_id}/{filler}/{pos}_neg_{neg}/benchmark.tsv"
#     params:
#         outdir = "out/meme/chip-diff_{meme_id}/{filler}/{pos}_neg_{neg}"
#     wildcard_constraints:
#         meme_id = "[a-zA-Z0-9-]+",
#         pos = "[a-zA-Z0-9_-]+",
#         neg = "[a-zA-Z0-9_-]+"
#     conda:
#         "../envs/meme.yaml"
#     shell:
#         """
#         meme-chip --oc {params.outdir} -neg {input.negfa} -db {input.meme} {input.posfa} &> {log}
#         """



rule meme_sea:
    """
    Doc:
        https://meme-suite.org/meme/doc/meme-chip.html?man_type=web
        http://meme-suite.org/doc/chip.html
    Warning:
        If the fasta input is a symlink, the tool will fail to appropriately symlink it again in its output directory...
    """
    input:
        fasta = "out/{filler}.fa",
        meme = lambda wildcards: mwconf['ids'][wildcards.meme_id]
    output:
        expand(
            "out/meme/sea_{{meme_id}}/{{filler}}/{file}",
            file=[
                'sea.html',
                'sea.tsv'
            ]
        )
    log:
               "out/meme/sea_{meme_id}/{filler}/log"
    benchmark:
               "out/meme/sea_{meme_id}/{filler}/benchmark.tsv"
    params:
        outdir="out/meme/sea_{meme_id}/{filler}"
    wildcard_constraints:
        meme_id = "[a-zA-Z0-9-]+"
    # sea install from conda is currently broken.
    #https://github.com/bioconda/bioconda-recipes/issues/22909
    # I have built it from source and installed it in ~/meme/bin/sea
    # Alternative would be to revert to conda version 5.0.2
    # conda:
        # "../envs/meme.yaml"
    shell:
        """
        ~/meme/bin/sea --oc {params.outdir} -m {input.meme} -p {input.fasta} &> {log}
        # sea --oc {params.outdir} -m {input.meme} -p {input.fasta} &> {log}
        """
