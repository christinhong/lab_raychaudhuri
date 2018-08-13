#' Get the variance explained by principal components
#' @param plist A list with class "prcomp" containing: sdev, rotation, x.
variance_explained <- function(plist) {
  rotation <- as.data.frame(plist$rotation)
  variance <- plist$sdev ^ 2  
  data.frame(
    Component = factor(1:ncol(rotation), levels = 1:ncol(rotation)),
    Variance = variance / sum(variance),
    CumulativeVariance = ( cumsum(variance) / sum(variance) )
  )
}

#' Plot the variance explained by the first few principal components.
#' 
#' @param plist A list with class "prcomp" containing: sdev, rotation, x,
#'  center, scale.
#' @param n The number of principal components to plot.
#' @return A ggplot2 plot.
plot_variance_explained <- function(plist, n = 10, cumulative = FALSE) {
  # Get the variance explained as a dataframe.
  dat <- head(variance_explained(plist), n)
  
  if (cumulative) {
    # Cumulative line plot.
    ggplot(data = dat) +
      geom_line(alpha = 0.8,
                aes(x = Component, y = CumulativeVariance, group = 1)) +
      geom_point(size = 3, aes(x = Component, y = CumulativeVariance)) +
      xlab("PC") +
      ylab("Cumulative Fraction of Variance Explained") +
      geom_abline(intercept = 0, slope = 100 * 1 / ncol(rotation),
                  alpha = 0.25) +
      scale_y_continuous(limits = c(0, max(dat$CumulativeVariance)))
  } else {
    # Bar plot.
    ggplot(data = dat) +
      geom_bar(alpha = 0.8, aes(x = Component, weight = Variance)) +
      ylab("Fraction of Variance Explained") +
      xlab("")
  }
}

#' Correlate principal components with factors in another data.frame.
#'
#' @param pca A "prcomp" object returned by prcomp().
#' @param df A matrix. Each column is tested for correlation with the PCs.
correlate_pcs = function(pca, df, npcs = 5, min.cor = 0.3) {
  pca.r = as.data.frame(pca$rotation)[ , 1:npcs]
  df = as.data.frame(df)
  # Make all of the columns numeric, so we can run cor().
  for (col in colnames(df)) {
    if (!class(df[ , col]) == "numeric") {
      df[ , col] = factorToNumeric(df[ , col])
    }
  }
  pca.cor = cor(cbind(pca.r, df), method = "spearman")
  result = list()
  for (col in colnames(df)) {
    # Exclude weak correlations.
    idx = abs(pca.cor[col, ]) > min.cor & abs(pca.cor[col, ]) < 1
    res = pca.cor[col, ][idx]
    res = res[!is.na(res)]
    if (length(res) > 1) {
      result[[col]] = res[order(abs(res), decreasing = TRUE)]
    }
  }
  result
}

#' List the genes with the greatest absolute values of loading factors.
#' @param pca A "prcomp" object returned by prcomp().
#' @param n_pcs Integer number of principal components to list.
#' @param n_items Integer number of items listed for each componenet.
#' @return A list of items with loading values.
loading_values = function(pca, n_pcs = 2, n_items = 30) {
  pc_scores = data.frame(pca$x)
  result = list()
  for (i in 1:n_pcs) {
    code = paste0("PC", i)
    sx = pc_scores[order(pc_scores[, code]), ]
    result[[code]] = c(
      head(sx[, code, drop = TRUE], n_items / 2),
      tail(sx[, code, drop = TRUE], n_items / 2)
    )
  }
  result
}

#' Return a numeric vector with the levels of a factor.
#'
#' The purpose of this function is to convert a non-numeric factor to a numeric
#' factor. This is useful when you want to compute correlation with non-numeric
#' vectors.
#'
#' @param xs A factor
#' @return A numeric vector.
factorToNumeric = function(xs) {
  uniqs = unique(xs)
  values = 1:length(uniqs)
  as.numeric(values[match(xs, uniqs)])
}
