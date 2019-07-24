def params_gemBS_bscall_e(wildcards):
    """
    Created:
        2018-04-04 19:19:02
    """
    r=wildcards['fa_genome_id']

    if r == 'GRCh38-BP-BS-protocol':
        e = 'HomoSapiens'
    
    return(e)

def input_gemBS_bscall_concatenate_bcf_list(wildcards):
    """
    Created:
        2018-04-04 19:30:15
                
    bcf="out/gemBS/bscall-concatenate_r-{fa_genome_id}{filler}_{chr}.bcf"

    """
    fa_genome_id = wildcards['fa_genome_id']
    filler = wildcards['filler']
    _p = wildcards['_p']

    if  fa_genome_id == 'GRCh38-BP-BS-protocol':
        paths = expand("out/gemBS/bscall"+_p+"_r-"+fa_genome_id+"/"+filler+"_{chr}.bcf", chr=HUMAN_MAIN_CHROMOSOMES_NCBI)

    return(paths)
    


