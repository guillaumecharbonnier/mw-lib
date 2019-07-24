####
### BIG TODO : parse config for bam requisites, and also the bai as it seems that when bai is missing I get this kind of error:
# ValueError: fetch called on bamfile without index
####

rule rgt_thor_tools_filter:
    """
    Created:
        2018-01-22 22:13:57
    Test:
        out/rgt/thor/test_thymus_H3K27ac/testName-pvalThreshold-5-diffpeaks.bed
    """
    input:
        zsh="opt/miniconda/envs/thor_tools/bin/zsh",
        filter="out/tar/xvzf/wget/http/www.regulatory-genomics.org/wp-content/uploads/2015/07/THOR-tools/filter-THOR.sh",
        diffpeaks="out/{filler}-diffpeaks.bed"
    output:
        diffpeaks="out/{filler}-pvalThreshold-{pvalThreshold}-diffpeaks.bed"
    shell:
        """
        {input.zsh} {input.filter} {input.diffpeaks} {wildcards.pvalThreshold} > {output.diffpeaks}
        """

rule rgt_thor_tools_split:
    """
    Created:
        2018-01-22 22:13:57
    Test:
        out/rgt/thor/test_thymus_H3K27ac/testName-pvalThreshold-5-diffpeaks-gain.bed
    """
    input:
        zsh="opt/miniconda/envs/thor_tools/bin/zsh",
        filter="out/tar/xvzf/wget/http/www.regulatory-genomics.org/wp-content/uploads/2015/07/THOR-tools/split-THOR.sh",
        diffpeaks="out/{filler}-diffpeaks.bed"
    output:
        diffpeaks_gain="out/{filler}-diffpeaks-gain.bed",
        diffpeaks_lose="out/{filler}-diffpeaks-lose.bed",
    shell:
        """
        {input.zsh} {input.filter} {input.diffpeaks}
        """


rule rgt_thor:
    """
    Created:
        2018-01-09 22:43:19
    Aim:
        Trying to produce normzalized bigwig for Salva and Vahid project.
    Test: 
        out/rgt/thor/test_thymus_H3K27ac/done
    Note:
        --name="THOR" so it does not add timestamp to output files which make 
        http://www.regulatory-genomics.org/wp-content/uploads/2015/07/THOR-tools.tar.gz
    """
    input:
        thor="opt/miniconda/envs/rgt/bin/rgt-THOR",
        config="src/rgt/thor/{design}.config"
    output:
        done="out/rgt/thor/{design}/done"
    params:
        outdir="out/rgt/thor/{design}"
    shell:
        """
        set +u; source opt/miniconda/bin/activate rgt; set -u
        
        rgt-THOR\
            --output-dir {params.outdir}\
            --name="THOR"\
            {input.config}
        """
        
rule rgt_thor_report:
    """
    Created:
        2018-01-09 22:43:19
    Aim:
        Trying to produce normzalized bigwig for Salva and Vahid project.
    Test: 
        out/rgt/thor_report/test_thymus_H3K27ac/done
    """
    input:
        thor="opt/miniconda/envs/rgt/bin/rgt-THOR",
        config="src/rgt/thor/{design}.config"
    output:
        done="out/rgt/thor_report/{design}/done"
    params:
        outdir="out/rgt/thor_report/{design}"
    shell:
        """
        set +u; source opt/miniconda/bin/activate rgt; set -u
        
        # rgt_thor_report failed because of this rerror:
        # IOError: [Errno 2] No such file or directory: '/cobelix/Charbonnier/rgtdata/fig/rgt_logo.gif'
        # I add this for the new test to see if it solves the issue.
        export RGTDATA={RGTDATA}

        rgt-THOR\
            --report\
            --output-dir {params.outdir}\
            {input.config}
        """

rule rgt_tmp:
    input:
        expand('out/rgt/thor_report_housekeeping_genes_deadzones/{config}/done', config=["H3K27ac_CD34_LC","H3K27ac_EC_LC","H3K27ac_SP4_LC","H3K27ac_SP8_LC"])


def input_rgt_thor_config_requirements(wildcards):
    """
    Created:
        2017-06-14 16:42:12 
    Aim:
        A function retrieving files in rgt thor config file.
    """
    id=wildcards['design']

    df = pandas.read_table("src/rgt/thor/"+id+".config", header=None, comment='#')
    print(df)
    # Column 2 should contain samples
    paths=df[0].tolist()
    print(paths)
    # 
    paths = [re.sub('.bam$','.bam.bai',x) for x in paths]
    return(paths)

rule rgt_thor_report_housekeeping_genes_deadzones:
    """
    Created:
        2018-01-09 22:43:19
    Aim:
        Trying to produce normalized bigwig for Salva and Vahid project.
    Test: 
        out/rgt/thor_report_housekeeping_genes_deadzones/test_thymus_H3K27ac/done
        out/rgt/thor_report_housekeeping_genes_deadzones/H3K27ac_CD34_EC/done
        out/rgt/thor_report_housekeeping_genes_deadzones/H3K27ac_CD34_EC_merged/done out/rgt/thor_report_housekeeping_genes_deadzones/H3K27ac_EC_LC_merged/done out/rgt/thor_report_housekeeping_genes_deadzones/H3K27ac_LC_SP4_merged/done out/rgt/thor_report_housekeeping_genes_deadzones/H3K27ac_LC_SP8_merged/done out/rgt/thor_report_housekeeping_genes_deadzones/H3K27ac_SP4_SP8_merged/done


    """
    input:
        thor="opt/miniconda/envs/rgt/bin/rgt-THOR",
        bed_hkg="out/bedtools/slop_g-hg38_b-1000/gtftk/5p_3p_coord/sed/add_chr/gtftk/select_by_key_key-transcript_id_file-with-values-human-hkg/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.bed",
        bed_blk="out/gunzip/to-stdout/wget/http/mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed",
        config="src/rgt/thor/{design}.config",
        config_requirements=input_rgt_thor_config_requirements
    output:
        done="out/rgt/thor_report_housekeeping_genes_deadzones/{design}/done"
    params:
        outdir="out/rgt/thor_report_housekeeping_genes_deadzones/{design}"
    shell:
        """
        set +u; source opt/miniconda/bin/activate rgt; set -u
        
        # rgt_thor_report failed because of this rerror:
        # IOError: [Errno 2] No such file or directory: '/cobelix/Charbonnier/rgtdata/fig/rgt_logo.gif'
        # I add this for the new test to see if it solves the issue.
        export RGTDATA={RGTDATA}

        rgt-THOR\
            --report\
            --output-dir {params.outdir}\
            --housekeeping-genes {input.bed_hkg}\
            --deadzones {input.bed_blk}\
            {input.config}
        
        touch {output.done}
        """

rule rgt_thor_report_housekeeping_genes:
    """
    Created:
        2018-01-09 22:43:19
    Aim:
        Trying to produce normzalized bigwig for Salva and Vahid project.
    Test: 
        out/rgt/thor_report_housekeeping_genes/test_thymus_H3K27ac/done
    """
    input:
        thor="opt/miniconda/envs/rgt/bin/rgt-THOR",
        bed_hkg="out/bedtools/slop_g-hg38_b-1000/gtftk/5p_3p_coord/sed/add_chr/gtftk/select_by_key_key-transcript_id_file-with-values-human-hkg/gunzip/to-stdout/wget/ftp/ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.bed",
        config="src/rgt/thor/{design}.config"
    output:
        done="out/rgt/thor_report_housekeeping_genes/{design}/done"
    params:
        outdir="out/rgt/thor_report_housekeeping_genes/{design}"
    shell:
        """
        set +u; source opt/miniconda/bin/activate rgt; set -u
        
        # rgt_thor_report failed because of this rerror:
        # IOError: [Errno 2] No such file or directory: '/cobelix/Charbonnier/rgtdata/fig/rgt_logo.gif'
        # I add this for the new test to see if it solves the issue.
        export RGTDATA={RGTDATA}

        rgt-THOR\
            --report\
            --output-dir {params.outdir}\
            --housekeeping-genes {input.bed_hkg}\
            {input.config}
        
        touch {output.done}
        """


