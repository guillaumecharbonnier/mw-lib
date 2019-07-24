rule test_gsea_with_known_genes_rnk:
	"""
	Created: 2016-07-08 11h05

	Testing how GSEA works with a totally not biological ranked list of genes.
	This rule take UCSC knownGenes ID and create a .rnk file for GSEA.	
	"""
	input:"annotation/input/mm9TSSRegions.bed"
	output:"out/gsea/test/knownGenesUniqTest.rnk"
	shell:"""
	cut -f4 annotation/input/mm9TSSRegions.bed | cut -f2 --delimiter='/' | sort -u | awk -F $'\t' 'BEGIN {OFS = FS}{print $0, NR}' > out/gsea/test/knownGenesUniqTest.rnk

	"""
