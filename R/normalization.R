#' Min-max normalize a matrix by columns
#'
#' This function performs min-max normalization for each column independently:
#' \deqn{(x - min(x)) / (max(x) - min(x))}
#'
#' Constant columns should be removed before calling this function.
#'
#' @param x A numeric matrix with samples/pixels as rows and features as columns.
#'
#' @return A numeric matrix of the same dimension as \code{x}, scaled to the
#' range \code{[0, 1]} for each column.
#'
#' @examples
#' mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
#' minmax_normalize(mat)
#'
#' @export
minmax_normalize <- function(x) {
  check_spectra_matrix(x)

  min_vals <- apply(x, 2, min, na.rm = TRUE)
  max_vals <- apply(x, 2, max, na.rm = TRUE)
  ranges <- max_vals - min_vals

  if (any(ranges == 0)) {
    stop(
      "Matrix contains constant columns. ",
      "Please remove them first using remove_constant_features()."
    )
  }

  x_scaled <- sweep(sweep(x, 2, min_vals, "-"), 2, ranges, "/")
  x_scaled <- as.matrix(x_scaled)

  rownames(x_scaled) <- rownames(x)
  colnames(x_scaled) <- colnames(x)

  x_scaled
}
