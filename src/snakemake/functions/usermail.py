def usermail(user):
    """
    Created:
        2019-01-21 11:09:48
    Aim:
        Change the destination of onsuccess email according to username.
    """
    d={}
    paths = glob.glob('../mw*/src/snakemake/tables/email*.tsv')
    #paths.extend(glob.glob('inp/*/src/snakemake/tables/email*.tsv'))
    #print(paths)
    for path in paths:
        #print('Looking for email in ' + path)
        fp = open(path)
        rdr = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
        for row in rdr:
            if row['user'] in d:
                print(row['user'] + ' from ' + path + ' is replacing previous value')
            d[row['user']] = row['email']

    return d.get(user)




#        fp = open(path)
#rdr = csv.DictReader(filter(lambda row: row[0]!='#', fp))
#print(rdr)
#for row in rdr:
#print(row)
#if row['alias'] in d:
#print(row['alias'] + ' from ' + path + ' is replacing previous value')
#d[row['alias']] = row['file']
#return d[alias_id]


