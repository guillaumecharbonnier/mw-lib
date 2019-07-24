#How to get a hash for a string:
#echo '--sortRegions no --interpolationMethod nearest --heatmapHeight 10 --heatmapWidth 1 --missingDataColor 1 --colorList '#313695,#B6DFEC,#FFF6B1,#F67C4A,#A70226' --whatToShow 'plot, heatmap and colorbar' --legendLocation upper-center--xAxisLabel ' ' --refPointLabel 0 --samplesLabel CD34 EC LC SP4 SP8' | sha1sum | cut -c -4

echo '--sortRegions no --heatmapHeight 10 --heatmapWidth 1 --colorList '#313695,#B6DFEC,#FFF6B1,#F67C4A,#A70226' --whatToShow 'plot, heatmap and colorbar' --legendLocation upper-center --xAxisLabel ' ' --refPointLabel 0 --samplesLabel HSC CD34 EC LC SP4 SP8' | sha1sum | cut -c -4

