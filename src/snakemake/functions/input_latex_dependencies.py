def latex_input_dependencies(wildcards):
    """
    Created:
        2018-10-14 03:52:13
    Modifed:
        2019-02-25 15:20:38 - Added pattern2 for amu_these. Ugly, could be generalized to any 'document_class' requiring cls.
    Aim:
        Parse content of a tex file to look for files in the same directory as the tex one that are dependencies. First need to retrieve the header and footer for invoice.
    """
    filler=wildcards['filler']
    tex="out/" + wildcards['filler'] + "/" + wildcards['filename'] + ".tex"
    print(tex)
    if os.path.isfile(tex):
        #pattern1 = re.compile(r"^[^%].*\\input.+{([^\}]+)}")
        pattern1 = re.compile(r"\\input{([^\}]+)}")
        pattern2 = re.compile(r"\\documentclass{amu_these}")
        paths = []
        
        with open (tex, "rt") as infile:
            for line in infile:
                m1 = pattern1.search(line)
                if m1:
                    path = "out/" + wildcards['filler'] + "/" + m1.group(1)
                    print(path)
                    paths.append(path)
                m2 = pattern2.search(line)
                if m2:
                    paths.append('out/' + wildcards['filler'] + "/amu_these.cls")
    else:
        paths=tex
    return(paths)

def latex_smi_dep(wildcards):
    filepaths = [] #Make sure returned paths is never empty.
    filepaths.append("out/" + wildcards['filler'] + "/" + wildcards['filename'] + ".tex")
    pattern = re.compile(r"^(?:(?!%).)*?\\def\\smi{([^\}]+)}")
    for filepath in filepaths:
        if filepath.endswith(".tex") or filepath.endswith(".cls"):
            print(filepath + ' is a tex file')
            if os.path.isfile(filepath):
                with open (filepath, "rt") as infile:
                    for line in infile:
                        match = pattern.search(line)
                        if match:
                            print(line)
                            filepath = match.group(1)
                            filepaths.append(filepath)
            else:
                eprint('The input tex file was not yet available for latex_smi_dep function to find input dependencies. You should try to run Snakemake again now.')

    return(filepaths)

def r_smi_dep(wildcards):
    filepaths = [] #Make sure returned paths is never empty.
    filepaths.append("out/" + wildcards['filler'] + ".R")
    pattern = re.compile(r"""smi *<- *["'](.+)["']""")
    #pattern = re.compile(r"""smi\s*<-\s*["'](.+)["']""")
    for filepath in filepaths:
        if filepath.endswith(".R"):
            print(filepath + ' is a R file')
            if os.path.isfile(filepath):
                with open (filepath, "rt") as infile:
                    for line in infile:
                        match = pattern.search(line)
                        if match:
                            print(line)
                            filepath = match.group(1)
                            filepaths.append(filepath)
            else:
                eprint('The input R file was not yet available for r_smi_dep function to find input dependencies. You should try to run Snakemake again now.')

    return(filepaths)

def r_smi_dir(dir):
    rmd_paths = glob.glob(dir + "/*.Rmd")
    dep_paths = []
    pattern = re.compile(r"""smi *<- *["'](.+)["']""")
    for rmd_path in rmd_paths:
        with open (rmd_path, "rt") as infile:
            for line in infile:
                match = pattern.search(line)
                if match:
                    print(line)
                    dep_paths.append(match.group(1))
    
    return(dep_paths)


#def parse_tex_file(path, pattern):
#    
#files_to_parse = ['toto', 'tata']
#for file_to_parse in files_to_parse:
#    print(file_to_parse)
#    files_to_parse.append('tutu')
#



def latex_includegraphics_dependencies_test(wildcards):
    """
    Created:
        2018-10-14 01:59:27
    Aim:
        A function parsing the content of a tex file and returning the paths of dependencies. It allows to write analysis directly from tex files without writing Snakefile.
    """
    filler=wildcards['filler']
    tex="out/" + wildcards['filler'] + "/" + wildcards['filename'] + ".tex"
    print(tex)
    if os.path.isfile(tex):
        # Meaning of the regex:
        # ^[^%].* means line that does not start with %, i.e. not commented in latex.
        # \\includegraphics.+{ means includegraphic call with options until start of file path.
        # (.*) is the file path extracted.
        # } is the end of the includegraphic call.
        # This regex implies avoiding additionnal {} in the file path.
        # pattern1 temporarily disabled for debugging - 2018-03-22 09:38:49
        #pattern1 = re.compile(r"IDontWantPattern1ToWorkCurrentlyButIDontWantToCommentCodeBelow")
        pattern1 = re.compile(r"^[^%].*\\include.+{([^\}]+)}")
        #pattern1 = re.compile(r"^[^%].*\\includegraphics.+{([^\}]+)}")
        # Note: pattern2 currently does not match if line start with \def. At least one character should be there first.
        pattern2 = re.compile(r"^[^%].*\\def\\snakemakeInput{([^\}]+)}")
        # Please help me finding the working regex for the two patterns at the same time:
        #pattern = re.compile(r"^[^%].*\\includegraphics.+{([^\}]+)}|^[^%].*\\def\\snakemakeInput{([^\}]+)}")
        # pattern3: 
        # Example:
        # \includepdf[pages=-]{out/gs/sDEVICE-pdfwrite_dCompatibilityLevel-1.5_dPDFSETTINGS-prepress/ln/alias/mendeley_path_to_bibtex_key/Goudarzi2016.pdf}
        # includepdf includegraphics
        
        #\def\snakemakeInput{doc/2017_01_06_gc_poster_hceres/screenshot/logos.pdf}
        paths = []
        
        with open (tex, "rt") as infile:
            for line in infile:
                m1 = pattern1.search(line)
                m2 = pattern2.search(line)
                #            m3 = pattern3.search(line)
                #m = pattern.search(line)
                if m1:
                    path = m1.group(1)
                    print(path)
                    if path == '\snakemakeInput':
                        print('snakemakeInput')
                    else:
                        paths.append(path)
                if m2:
                    path = m2.group(1)
                    print('snakemakeInput path: ' + path)
                    paths.append(path)
                #            if m3:
                #                print('m3 found')
                #                path = m3.group(1)
                #                paths.append(path)
    else:
        paths=tex
    return(paths)


def latex_includegraphics_dependencies(wildcards):
    """
    Created:
        2017-01-31 15:18:34
    Aim:
        A function parsing the content of a tex file and returning the paths of dependencies. It allows to write analysis directly from tex files without writing Snakefile.
    """
    id=wildcards['id']
    tex="src/tex/"+id+".tex"
    
    # Meaning of the regex:
    # ^[^%].* means line that does not start with %, i.e. not commented in latex.
    # \\includegraphics.+{ means includegraphic call with options until start of file path.
    # (.*) is the file path extracted.
    # } is the end of the includegraphic call.
    # This regex implies avoiding additionnal {} in the file path.
    # pattern1 temporarily disabled for debugging - 2018-03-22 09:38:49
    #pattern1 = re.compile(r"IDontWantPattern1ToWorkCurrentlyButIDontWantToCommentCodeBelow")
    pattern1 = re.compile(r"^[^%].*\\include.+{([^\}]+)}")
    #pattern1 = re.compile(r"^[^%].*\\includegraphics.+{([^\}]+)}")
    pattern2 = re.compile(r"^[^%].*\\def\\snakemakeInput{([^\}]+)}")
    # Please help me finding the working regex for the two patterns at the same time:
    #pattern = re.compile(r"^[^%].*\\includegraphics.+{([^\}]+)}|^[^%].*\\def\\snakemakeInput{([^\}]+)}")
    # pattern3: 
    # Example:
    # \includepdf[pages=-]{out/gs/sDEVICE-pdfwrite_dCompatibilityLevel-1.5_dPDFSETTINGS-prepress/ln/alias/mendeley_path_to_bibtex_key/Goudarzi2016.pdf}
    # includepdf includegraphics

    #\def\snakemakeInput{doc/2017_01_06_gc_poster_hceres/screenshot/logos.pdf}
    paths = []

    with open (tex, "rt") as infile:
        for line in infile:
            m1 = pattern1.search(line)
            m2 = pattern2.search(line)
#            m3 = pattern3.search(line)
            #m = pattern.search(line)
            if m1:
                path = m1.group(1)
                print(path)
                if path == '\snakemakeInput':
                    print('snakemakeInput')
                else:
                    paths.append(path)
            if m2:
                path = m2.group(1)
                paths.append(path)
#            if m3:
#                print('m3 found')
#                path = m3.group(1)
#                paths.append(path)
#
    return(paths)

