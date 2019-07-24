def input_yaml_greatr(wildcards):
    """
    Created:
        2019-01-24 14:12:39
    """
    # Getting variables from Snakemake wildcards.
    bed_list_id=wildcards['bed_list_id']
    d={}
    #paths = glob.glob('../mw*/src/snakemake/tables/yaml_greatr*.tsv')
    path = glob.glob('../mw*/src/greatr/' + bed_list_id + '.yaml')
    if len(path) == 1:
        #return(path)
        print('Loading greatr configuration file: ' + path[0])
    else:
        print('More than one yaml exists for this bed list:')
        print(path)
        print('The first one will be used.')
    return(path[0])
    #for path in paths:
    #    print('Looking for yaml greatr id in ' + path)
    #    fp = open(path)
    #    rdr = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
    #    for row in rdr:
    #        if row['bed_list_id'] in d:
    #            print(row['bed_list_id'] + ' from ' + path + ' is replacing previous value')
    #        d[row['bed_list_id']] = row['yaml_path']
    #return d[bed_list_id]

