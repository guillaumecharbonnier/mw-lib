# Note: color_palettes is defined in manual_palettes.R (sourced earlier).
# Do NOT redefine it here as it would override the comprehensive palette from manual_palettes.R.

plotPcaIndFacetsForEachCov <- function(
                                               d = d_pca,
                                               d_cov = d_covariates,
                                               n_row = 1,
                                               n_col = ncol(d_covariates),
                                               axes_to_plot = c(1, 2)) {
  # Align d_cov to active PCA individuals (handles spectator designs where
  # d_cov may contain supplementary samples not used in prcomp)
  if (!is.null(rownames(d$x)) && !is.null(rownames(d_cov))) {
    d_cov <- d_cov[rownames(d$x), , drop = FALSE]
  }
  plots <- list()
  for (covariate_name in names(d_cov)) {
    plot <- fviz_pca_ind(
      d,
      axes = axes_to_plot,
      mean.point = FALSE,
      label = "none",
      habillage = d_cov[[covariate_name]],
      legend.title = covariate_name
    )
    plot <- plot + ggtitle(NULL)
    # Apply custom colors from manual_palettes.R if available for this covariate
    if (covariate_name %in% names(color_palettes)) {
      plot <- plot + scale_color_manual(values = color_palettes[[covariate_name]])
    }
    plots[[covariate_name]] <- plot
  }

  ggarrange(
    plotlist = plots,
    ncol = n_col,
    nrow = n_row,
    align = "hv",
    common.legend = FALSE
  )
}

plotUmapIndFacetsForEachCov <- function(
  d = d_umap,
  d_cov = d_covariates,
  n_row = 1,
  labels = NULL,
  n_col = ncol(d_covariates)
) {
  d_umap <- data.table(
    data.table(d$layout),
    as.data.table(d_cov)
  )

  pl <- list()
  plots <- list()
  for (covariate_name in names(d_cov)) {
    plot <- ggplot(
      d_umap,
      aes(
      x = V1,
      y = V2,
      color = .data[[covariate_name]],
      shape = .data[[covariate_name]]
      )
    )
    plot <- plot + scale_shape_manual(values=1:nlevels(factor(d_cov[[covariate_name]])) + 1)
    if (!is.null(labels)) {
      plot <- plot + geom_text(aes(label = labels), size = 3)
    } else {
      plot <- plot + geom_point()
    }
    # Add custom colors for covariates with defined palette
    if (covariate_name %in% names(color_palettes)) {
      plot <- plot + scale_color_manual(values = color_palettes[[covariate_name]])
    }
    # Add minimal theme
    plot <- plot + theme_minimal()

    # Add UMAP1 and UMAP2 as axis labels
    plot <- plot + labs(x = "UMAP1", y = "UMAP2")

    #plot <- plot + ggtitle(NULL)
    plots[[covariate_name]] <- plot
  }

  ggarrange(
    plotlist = plots,
    ncol = n_col,
    nrow = n_row,
    align = "hv",
    common.legend = FALSE
  )
}
