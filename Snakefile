# Default metaworkflow Snakefile

import os
#Needed to have access to {WDIR} inside rules
WDIR= os.getcwd()
import glob

# Including other python imports.
paths = glob.glob("../mw*/src/snakemake/imports.py")

# Including global constraints on wildcards
paths.extend(glob.glob("../mw*/src/snakemake/wildcard_constraints.smk"))

# Including the list of ids and parameters
paths.extend(glob.glob("../mw*/src/snakemake/variables.py"))
paths.extend(glob.glob("../mw*/src/snakemake/variables.snake"))

# Including functions to map files to ids.
paths.extend(glob.glob("../mw*/src/snakemake/functions/*.py"))
paths.extend(glob.glob("../mw*/src/snakemake/functions/*.snake"))

# Including rules
paths.extend(glob.glob('../mw*/src/snakemake/rules/*.smk'))
paths.extend(glob.glob('../mw*/src/snakemake/rules/*.rules'))

for path in paths:
    include: path
    #eprint("Loaded: " + path)

# ids may already be filled with some function like in mw-sst plugin
if 'ids' not in config:
    config['ids']={}

# Loading config dicts
paths = glob.glob('../mw*/src/snakemake/tables/*ids.tsv')

for path in paths:
    #print('Looking for alias in ' + path)
    fp = open(path)
    rdr = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
    for row in rdr:
        if row['id'] in config["ids"]:
            eprint(row['id'] + ' from ' + path + ' is replacing previous value')
        config["ids"][row['id']] = row['path']

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input:
        #Maybe add here the compilated doc for metaworkflow?

MAILAPP="mail" #or "mutt" for laptop if configured with password. Test that later.
USER=getpass.getuser()
USERMAIL = usermail(USER)

paths_onstart = glob.glob('../mw*/src/snakemake/hooks/onstart/*')
onstart:
    shell("for path in {paths_onstart}; do . $path; done")

paths_onsuccess = glob.glob('../mw*/src/snakemake/hooks/onsuccess/*')
onsuccess:
    shell("for path in {paths_onsuccess}; do . $path; done")
    shell("{MAILAPP} -s 'workflow finished' {USERMAIL} < {log}")
# mutt -s "Test from mutt" {USERMAIL} < {log}

paths_onerrror = glob.glob('../mw*/src/snakemake/hooks/onerror/*')
onerror:
    shell("for path in {paths_onerror}; do . $path; done")
    shell("{MAILAPP} -s 'an error occurred' {USERMAIL} < {log}")

#Uncomment these for debugging:
#print(dir()) #will give you the list of in scope variables:
#print(globals()) #will give you a dictionary of global variables
#print(locals())

