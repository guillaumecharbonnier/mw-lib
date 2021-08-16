plotPcaIndFacetsForEachCov <- function(
                                               d = d_pca,
                                               d_cov = d_covariates,
                                               n_row = 1,
                                               n_col = ncol(d_covariates),
                                               axes_to_plot = c(1, 2)) {
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