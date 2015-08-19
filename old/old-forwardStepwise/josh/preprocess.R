# Functions for preprocessing design matrices

#' Encode factors as groups of dummy variables.
#'
#' @param x Design matrix
#' @param index Group membership indicator of length p
#' @param cols Indicator of columns to be converted to groups
#' @param center Centers group of dummy variables, default TRUE
#' @return Expanded design matrix with dummy variables
factor_matrix <- function(x, index, cols, center = FALSE, nobs.thresh = NULL, ...) UseMethod("factor_matrix")

factor_matrix.default <- function(x, index, cols, center = FALSE, nobs.thresh = NULL, ...) {

  if (missing(index)) index <- 1:ncol(x)
  if (missing(cols)) cols <- 1:ncol(x)

  ldata <- lapply(1:ncol(x), function(j) {
    if (j %in% cols) {
      # Categorical, create dummy variables using model.matrix

      if (((is.null(nobs.thresh)) || (min(table(x[, j])) >= nobs.thresh)) && (length(unique(x[, j])) > 1)) {
        # Omit categories with rare levels if nobs.thresh is defined.

        submat <- model.matrix(~ factor(x[, j]) - 1)
        if (center) submat <- submat - mean(submat)
        colnames(submat) <- paste(colnames(x)[j], 1:ncol(submat), sep=":")

        x_out <- submat
        index_out <- rep(index[j], ncol(submat))

      } else {
        x_out <- as.vector(rep(NA, nrow(x)))
        index_out <- 0
      }
    } else {
      # Not categorical. Keep without changing.
      x_out <- x[, j]
      index_out <- index[j]

    }
    if (length(index_out) == 1) return(matrix(c(index_out, x_out), ncol=1))
    return(rbind(index_out, x_out))
  })

  out <- do.call(cbind, ldata)
  out <- out[, which(!is.na(out[nrow(out), ]))]
  xout <- out[-1, ]
  if (!is.null(colnames(x))) attr(xout, "varnames") <- colnames(x)

  invisible(list(x=xout, index=out[1, ]))
}


#' Scale design matrix by groups
#'
#' @param x Design matrix
#' @param index Group membership indicator of length p
#' @param center Center groups, default TRUE
#' @param scale Scale groups by Frobeniusm norm, default TRUE
#' @return Scaled design matrix
scale_groups <- function(x, index, center = TRUE, scale = TRUE, ...) UseMethod("scale_groups")

scale_groups.default <- function(x, index, center = TRUE, scale = TRUE, ...) {

  for (j in unique(index)) {
    inds <- which(index == j)
    if (center) x[, inds] <- x[, inds] - mean(x[, inds])
    if (scale) x[, inds] <- x[, inds] / sqrt(sum(x[, inds]^2))

  }
  return(x)
}

#' Shift index labels to 1:G, where G is the number of groups
#'
#' @param index Current index indicator
#' @return List of length two, first the new index and second the labels
shift_index <- function(index) {
    lengths <- rle(index)$lengths
    labels <- 1:length(lengths)
    index <- rep(labels, lengths)
    return(list(index = index, labels = labels))
}
