rule deepTools_plotProfile_extra:
    """
    Created:
        2017-10-23 17:12:24
    Test:
        out/texlive/pdfcrop/deepTools/plotProfile_--numPlotsPerRow_1_--refPointLabel_0_--samplesLabel_thymus-stages/deepTools/computeMatrix_reference-point_--referencePoint_center_-b_5000_-a_5000_-bs_200_--sortRegions_keep_-R_hg38-dClust-rowFeature-no-rmsk-mxy-no-donor-effect-distal_-S_hg38-H3K27ac-thymus-merged-wiq.pdf
    """
    input:
        matrix="out/{filler}.txt.gz",
    output:
        pdf="out/{tool}{extra}/{filler}.pdf"
    params:
        extra = params_extra
    threads:
        1
    wildcard_constraints:
        tool = "deepTools/plotProfile"
    conda:
        "../envs/deeptools.yaml"
    priority:
        2
    shell:
        "plotProfile --matrixFile {input.matrix} --outFileName {output.pdf} {params.extra}"

#####################
# ONLY LEGACY BELOW #
#####################

rule deepTools_plotProfile_perGroup:
    """
    Created:
        2017-10-23 17:12:24
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotProfile="opt/miniconda/envs/deeptools/bin/plotProfile"
    output:
        pdf="out/deepTools/plotProfile_perGroup/{filler}.pdf",
    threads:
        1
    priority:
        2
    shell:
        """
        {input.plotProfile} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --perGroup
        """

rule deepTools_plotProfile_kmeans:
    """
    Created:
        2017-10-23 17:12:24
    Test:
    """
    input:
        matrix="out/{filler}.txt.gz",
        plotProfile="opt/miniconda/envs/deeptools/bin/plotProfile"
    output:
        pdf="out/deepTools/plotProfile_kmeans-{kmeans}/{filler}.pdf"
    threads:
        1
    priority:
        2
    shell:
        """
        {input.plotProfile} \
            --matrixFile {input.matrix} \
            --outFileName {output.pdf} \
            --kmeans {wildcards.kmeans}
        """

rule deepTools_plotProfile_legacy:
	""" 
	Created: 2016-03-30
	
	Usage:
		expand("result/deepTools/plotProfile/scaleRegions/{signal}/{index}/{exp}/{feature}/{id}.pdf", signal="bamCompare", index="mm10", exp="merge_run113_run125_run126", feature="knownGene", id=["MNS-SC-WT_Log2ratioOver_MNS-SC-KO", "PSK-SC-WT_Log2ratioOver_PSK-SC-KO"])

	"""
	input:	matrix="out/deepTools/computeMatrix/{matrix_type}/{signal}/{index}/{exp}/{feature}/{lengthAround}bp/{id}.txt.gz",\
		plotProfile="opt/miniconda/envs/py27/bin/plotProfile"
	output: pdf="result/deepTools/plotProfile/{matrix_type}/{signal}/{index}/{exp}/{feature}/{lengthAround}bp/{id}.pdf"
	threads: 1
	shell:"""
	{input.plotProfile} \
		--matrixFile {input.matrix}\
		--outFileName {output.pdf}
	"""	

