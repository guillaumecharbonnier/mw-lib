# This function solves limitations in default aheatmap() palette for annCol.
# Provide as input the data.frame you provide for aheatmap annCol arg
# and it return the list object you should provide to aheatmap annColors.
#annColToAnnColors <- function(ann_col){
#    ann_colors <- list()
#    pal <- brewer.pal(n = 12, name = "Set3")
#
#    n_quali_vars <- length(unlist(apply(ann_col,2,unique)))
#    
#    pal <- createPalette(n_quali_vars, brewer.pal(n = 8, name = "Set2"))
#    names(pal) <- NULL
#    #pal <- colorRampPalette(pal)(n_quali_vars)
#
#    idx_start <- 1
#    for (factor_col_name in colnames(ann_col)){
#        idx_end <- idx_start + length(levels(ann_col[,factor_col_name]))
#        ann_colors[[factor_col_name]] <- pal[idx_start:idx_end]
#        idx_start <- idx_end 
#    }
#    return(ann_colors)
#}

annColToAnnColors <- function(ann_col){
    ann_colors <- list()
    seed_pal <- brewer.pal(n = 9, name = "Set1")

    pal_start <- 0
    for (factor_col_name in colnames(ann_col)){
        pal <- createPalette(length(levels(ann_col[,factor_col_name])),
                             seed_pal[pal_start+1])
        names(pal) <- NULL
        ann_colors[[factor_col_name]] <- pal
        pal_start <- (pal_start + 1) %% 9
    }
    return(ann_colors)
}

