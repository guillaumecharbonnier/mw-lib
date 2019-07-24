def params_deepTools_plotHeatmap_colorList(wildcards):
    """
    Created:
        2017-05-07 18:26:4
    Aim: 
        Replace palette function with new rules for paths.
    """
    tmp = str(wildcards['deepTools_plotHeatmap_colorList_id']);
    
    #v2.2 style
    #colorList = {
    #    'whiteRed':                 "'#FFFFFF' '#EF0214'",
    #    'yellowOrangeRed':          "'#FFF6B1' '#F67C4A' '#A70226'",
    #    'blueCyanYellowOrangeRed':  "'#313695' '#B6DFEC' '#FFF6B1' '#F67C4A' '#A70226'",
    #    'blueYellowRed':            "blue yellow red",
    #    'blueLightyellowRed':       "blue lightyellow red",
    #    'cornflowerblueLightyellowIndianred': "cornflowerblue lightyellow indianred",
    #    'whiteIndianred':           "white indianred",
    #    'blueYRed':                 "blue y red"};
    
    #v2.4 style
    colorList = {
        'whiteRed':                 "'#FFFFFF,#EF0214'",
        'yellowOrangeRed':          "'#FFF6B1,#F67C4A,#A70226'",
        'blueCyanYellowOrangeRed':  "'#313695,#B6DFEC,#FFF6B1,#F67C4A,#A70226'",
        'blueYellowRed':            "blue,yellow,red",
        'blueLightyellowRed':       "blue,lightyellow,red",
        'cornflowerblueLightyellowIndianred': "cornflowerblue,lightyellow,indianred",
        'whiteIndianred':           "white,indianred",
        'blueYRed':                 "blue,y,red"};

    return "--colorList " + colorList[tmp]

def params_deepTools_plotHeatmap_whatToShow(wildcards):
    """
    Created:
        2017-05-07 18:36:17
    Aim:
        what to show take argument with spaces which is not very good with paths.
    """
    tmp = str(wildcards['deepTools_plotHeatmap_whatToShow_id']);

    whatToShow = {
        'phc': "'plot, heatmap and colorbar'",
        'pc':  "'plot and heatmap'",
        'h':   "'heatmap only'",
        'hc':  "'heatmap and colorbar'"
        }

    return "--whatToShow " + whatToShow[tmp]

def params_deepTools_plotHeatmap_xAxisLabel(wildcards):
    """
    Created:
        2017-05-07 18:51:43
    Aim:
        what to show take argument with spaces which is not very good with paths.
        xAxisLabel-tss/xAxisLabel-peak-center
    """
    tmp = str(wildcards['deepTools_plotHeatmap_xAxisLabel_id']);

    xAxisLabel = {
        'tss': "'TSS'",
        'gene-distance-bp': "'gene distance'",
        'nuc-center':  "'nuc. center'",
        'peak-center':  "'peak center'",
        'distance-to-nucleosome-center-bp':  "'nuc. distance to nucleosome center'"
        }

    return "--xAxisLabel " + xAxisLabel[tmp]

def params_deepTools_plotHeatmap_refPointLabel(wildcards):
    """
    Created:
        2017-05-07 18:58:06
    Aim:
        what to show take argument with spaces which is not very good with paths.
    """
    tmp = str(wildcards['deepTools_plotHeatmap_refPointLabel_id']);

    refPointLabel = {
        'TSS': "'TSS'",
        '0':  "'0'"
        }

    return "--refPointLabel " + refPointLabel[tmp]

def params_deepTools_plotHeatmap_startLabel(wildcards):
    """
    Created:
        2017-08-02 13:14:14
    Aim:
        what to show take argument with spaces which is not very good with paths.
    """
    tmp = str(wildcards['deepTools_plotHeatmap_startLabel_id']);

    startLabel = {
        'start': "'start'",
        '0':  "'0'"
        }

    return "--startLabel " + startLabel[tmp]

def params_deepTools_plotHeatmap_endLabel(wildcards):
    """
    Created:
        2017-08-02 13:15:45
    Aim:
        what to show take argument with spaces which is not very good with paths.
    """
    tmp = str(wildcards['deepTools_plotHeatmap_endLabel_id']);

    endLabel = {
        'end': "'end'",
        '0':  "'0'"
        }

    return "--endLabel " + endLabel[tmp] 

def params_deepTools_plotHeatmap_palette(wildcards):
    """
    Created: 2017-01-09 10h160 - A fonction to define palettes used as argument by deepTools_plotHeatmap.
    Should replace any specific palette function defined below for testing.
    """
    tmp = str(wildcards['palette']);
    colors = {
        'bwr': '--colorMap bwr',
        'customHexCol_whiteRed': "--colorList '#FFFFFF' '#EF0214'",
        # The way colorList work has changed in recent version of DeepTools:
        # 2.2.4 style:
        'customHexCol_YellowOrangeRed': "--colorList '#FFF6B1' '#F67C4A' '#A70226'",
        # 2.4.0 style:
        #'customHexCol_YellowOrangeRed': "--colorList '#FFF6B1,#F67C4A,#A70226'",
        # 2.2.4 style:
        #'customHexCol_blueCyanYellowOrangeRed': "--colorList '#313695' '#B6DFEC' '#FFF6B1' '#F67C4A' '#A70226'",
        # 2.4.0 style:
        'customHexCol_blueCyanYellowOrangeRed': "--colorList '#313695,#B6DFEC,#FFF6B1,#F67C4A,#A70226'",
        'blueYellowRed': '--colorList blue yellow red',
        'blueLightyellowRed': '--colorList blue lightyellow red',
        'cornflowerblueLightyellowIndianred': '--colorList cornflowerblue lightyellow indianred',
        'whiteIndianred': '--colorList white indianred',
        'RdYlBu': '--colorMap RdYlBu',
        'blueYRed': '--colorList blue y red'};
    return colors[tmp] 

#Probably obsolete with new way of writing rules.
#def deepTools_plotHeatmap_clustering(wildcards):
#    """
#    Created: 2017-01-09 11h16 - A fonction to define clustering approach.
#    """
#    tmp = str(wildcards['clustering']);
#    choices = {
#        'noClustering': '',
#        'kmeans': "CHANGETHAT",
#        'ADDOTHERCLUSTERINGAVAILABLEHERE': ''};
#    return choices[tmp]


