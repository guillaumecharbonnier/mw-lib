def params_extra(wildcards):
    """
    Created:
        2018-10-22 15:41:14
    Aim:
        Take the extra wildcard and return it formatted as a string for command line.
        It converts '_' to ' ' the one listed in 'src/snakemake/tables/protectedunderscores.tsv'. These string could contain argument names or values which naturally contain '_'.
        It also convert argument ids to their values accepted by tools as described in 'src/snakemake/tables/argids.tsv'.
    """
    extra=wildcards['extra']

    paths = glob.glob('../mw*/src/snakemake/tables/extra/*protectedunderscores.tsv')
    for path in paths:
        #print(path)
        for row in csv.DictReader(open(path), delimiter='\t', quoting=csv.QUOTE_NONE):
            if row['tool'] == wildcards['tool']:
                strings_to_protect = row['protectedunderscore'].split(",")
                print(strings_to_protect)
                for string_to_protect in strings_to_protect:
                    protected_string = re.sub('_','protectedunderscore',string_to_protect)
                    extra = re.sub(string_to_protect,protected_string,extra)

    extra=re.sub('_',' ',extra)
    extra=re.sub('protectedunderscore','_',extra)

    paths = glob.glob('../mw*/src/snakemake/tables/extra/argids.tsv')
    for path in paths:
        for row in csv.DictReader(open(path),delimiter='\t', quoting=csv.QUOTE_NONE):
            if row['tool'] == wildcards['tool']:
                extra = re.sub(row['arg'] + ' ' + row['id'], row['arg'] + ' ' + row['value'], extra)

    paths = glob.glob('../mw*/src/snakemake/tables/extra/extraids.tsv')
    for path in paths:
        for row in csv.DictReader(open(path),delimiter='\t', quoting=csv.QUOTE_NONE):
            if row['tool'] == wildcards['tool']:
                extra = re.sub(row['id'], row['extra'], extra)

    return(extra)
