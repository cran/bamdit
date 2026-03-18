#' Basic function to plot the data of meta-analysis of diagnostic test
#'
#' This function plots the true positive rates vs the false positive rates of each study included
#' in the meta-analysis. Study results are displayed by circles, the diameter of each circle is proportional
#' to the sample size of the study (or table). If subgroups are displayed each group is represented by
#' different colours. This function use the package \emph{ggplot2}.
#'
#'
#' @param data Either a data frame with at least 4 columns containing the true positives (tp),
#' number of patients with disease (n1), false positives (fp), number of patients without
#' disease (n2), or for two.by.two = TRUE a data frame where each line contains the
#' diagnostic results as a two by two table, where the column names are:
#' TP, FP, TN, FN.
#' @param two.by.two If TRUE indicates that the diagnostic results are given as: TP, FP, TN, FN.
#' @param group a variable name indicating a group factor
#' @param x.lo lower limit of the x-axis
#' @param x.up upper limit of the x-axis
#' @param y.lo lower limit of the y-axis
#' @param y.up upper limit of the y-axis
#' @param alpha.p transparency of the points
#' @param max.size scale parameter of the maximum size
#' @param points.color color of the points when no grouping variable is provided (default = "royalblue")
#' @param group.colors optional vector of colors for groups. If fewer colors than groups are provided,
#' the function will automatically generate additional colors using interpolation.
#'
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#'
#' data(ct)
#' ct$design <- with(ct, factor(design,
#'              labels = c("Prospective", "Retrospective")))
#'
#' plotdata(ct,              # Data frame
#'         group = "design", # Groupping variable
#'         y.lo = 0.75,      # Lower limit of y-axis
#'         x.up = 0.75,      # Upper limit of x-axis
#'         alpha.p = 0.5,    # Transparency of the balls
#'         max.size = 5)     # Scale the circles
#'}
#'
#'
#' @import ggplot2
#' @export
plotdata <- function(data,
                     two.by.two = FALSE,
                     group = NULL,
                     x.lo = 0, x.up = 1,
                     y.lo = 0, y.up = 1,
                     alpha.p = 0.7,
                     max.size = 15,
                     points.color = "royalblue",
                     group.colors = NULL)
{ gr=NULL
  if(two.by.two == FALSE)
  {
    tp <- data[,1]
    n1 <- data[,2]
    fp <- data[,3]
    n2 <- data[,4]
  } else
  {
    tp <- data$TP
    fp <- data$FP
    fn <- data$FN
    tn <- data$TN
    n1 <- tp + fn
    n2 <- fp + tn
  }

  # Data errors
  if(any(tp > n1) || any(fp > n2)) stop("the data is inconsistent")

  if(!missing(data)){
    tpr <- tp / n1
    fpr <- fp / n2
    n <- n1 + n2
  } else {
    stop("NAs are not alow in this plot function")
  }

  # ---- Build plotting data ----
  if(is.null(group)){
    dat.plot <- data.frame(tpr, fpr, n)
  } else {
    dat.plot <- data.frame(tpr, fpr, n, gr = data[[group]])
  }

  # ---- Plot with grouping ----
  if(!is.null(group)){

    p <- ggplot(dat.plot, aes(x = fpr, y = tpr, size = n, fill = gr)) +
      scale_x_continuous(name = "FPR (1 - Specificity)", limits = c(x.lo, x.up)) +
      scale_y_continuous(name = "TPR (Sensitivity)", limits = c(y.lo, y.up)) +
      geom_point(shape = 21, alpha = alpha.p) +
      scale_size_area(max_size = max.size)

    # ---- Apply custom group colors if provided ----
    if(!is.null(group.colors)){

      n_groups <- length(unique(dat.plot$gr))

      # If fewer colors than groups, interpolate additional colors
      if(length(group.colors) < n_groups){
        group.colors <- grDevices::colorRampPalette(group.colors)(n_groups)
      } else {
        group.colors <- group.colors[seq_len(n_groups)]
      }

      p <- p + scale_fill_manual(values = group.colors)
    }

    return(p)

  } else {

    # ---- Plot without grouping ----
    ggplot(dat.plot, aes(x = fpr, y = tpr, size = n)) +
      scale_x_continuous(name = "FPR (1 - Specificity)", limits = c(x.lo, x.up)) +
      scale_y_continuous(name = "TPR (Sensitivity)", limits = c(y.lo, y.up)) +
      geom_point(shape = 21, fill = points.color, alpha = alpha.p) +
      scale_size_area(max_size = max.size)
  }
}
