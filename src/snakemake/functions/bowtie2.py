def input_bowtie2_index_parts(wildcards):
    """
    Created:
        2018-10-22 15:41:14
    Aim:
        Take the extra wildcard and return it formatted as a string for command line.
        It converts '_' to ' ' the one listed in 'src/snakemake/tables/protectedunderscores.tsv'. These string could contain argument names or values which naturally contain '_'.
        It also convert argument ids to their values accepted by tools as described in 'src/snakemake/tables/argids.tsv'.
    """
    for row in csv.DictReader(open('../mw-lib/src/snakemake/tables/bowtie2_index_parts.tsv'), delimiter='\t', quoting=csv.QUOTE_NONE):
        if row['assembly'] == wildcards['index']:
            paths = eval(row['bowtie2_index_parts'])
    return(paths)

def params_bowtie2_index_base_path(wildcards):
    basepath = input_bowtie2_index_parts(wildcards)
    basepath = os.path.commonprefix(basepath)[:-1] # We do not want the last '.'
    return(basepath)


