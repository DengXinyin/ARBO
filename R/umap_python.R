#' Run UMAP using a Python backend via reticulate
#'
#' This function calls Python \pkg{umap-learn} through \pkg{reticulate} and
#' returns the embedding as an R data.frame.
#'
#' A working Python environment with \pkg{umap-learn} and \pkg{numpy} installed
#' is required. This function was tested with Python 3.9, \pkg{numpy} 1.24.4,
#' and \pkg{umap-learn} 0.5.7 in a conda environment, although other compatible
#' versions may also work.
#'
#' @param x A numeric matrix with samples/pixels as rows and features as columns.
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#' @param metric Distance metric used by UMAP. Default is \code{"cosine"}.
#' @param n_neighbors Integer; number of neighbors used by UMAP. Default is \code{10L}.
#' @param min_dist Numeric; UMAP \code{min_dist} parameter. Default is \code{0.05}.
#' @param n_components Integer; number of embedding dimensions. Default is \code{2L}.
#' @param random_state Integer; random seed for UMAP. Default is \code{2025L}.
#' @param n_jobs Integer; number of parallel jobs used by UMAP. Default is \code{1L}.
#' @param verbose Logical; whether to print UMAP progress. Default is \code{TRUE}.
#'
#' @return A data.frame containing UMAP coordinates. Column names are
#' \code{UMAP1}, \code{UMAP2}, ... depending on \code{n_components}.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(runif(100), nrow = 10)
#' emb <- run_umap_py(mat, python_path = "/path/to/python")
#' head(emb)
#' }
#'
#' @export
run_umap_py <- function(
    x,
    python_path = NULL,
    metric = "cosine",
    n_neighbors = 10L,
    min_dist = 0.05,
    n_components = 2L,
    random_state = 2025L,
    n_jobs = 1L,
    verbose = TRUE
) {
  check_spectra_matrix(x)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }

  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  umap <- reticulate::import("umap", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

  np$random$seed(as.integer(random_state))

  umap_result <- umap$UMAP(
    metric = metric,
    n_neighbors = as.integer(n_neighbors),
    min_dist = min_dist,
    n_components = as.integer(n_components),
    random_state = as.integer(random_state),
    n_jobs = as.integer(n_jobs),
    verbose = verbose
  )$fit_transform(x)

  umap_matrix <- reticulate::py_to_r(
    np$array(umap_result, dtype = np$float64)
  )

  umap_df <- as.data.frame(umap_matrix)
  colnames(umap_df) <- paste0("UMAP", seq_len(n_components))
  rownames(umap_df) <- rownames(x)

  umap_df
}


#' Check Python UMAP availability
#'
#' This helper function checks whether the required Python packages
#' \pkg{umap-learn} and \pkg{numpy} are available in the current or specified
#' Python environment.
#'
#' The Python backend was tested with Python 3.9, \pkg{numpy} 1.24.4, and
#' \pkg{umap-learn} 0.5.7 in a conda environment, although other compatible
#' versions may also work.
#'
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#'
#' @return Invisibly returns \code{TRUE} if Python UMAP dependencies are available.
#'
#' @examples
#' \dontrun{
#' check_umap_python_env("/path/to/python")
#' }
#'
#' @export
check_umap_python_env <- function(python_path = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }

  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }

  py_config <- reticulate::py_config()

  if (!reticulate::py_module_available("umap")) {
    stop("Python module 'umap' is not available in the selected environment.")
  }

  if (!reticulate::py_module_available("numpy")) {
    stop("Python module 'numpy' is not available in the selected environment.")
  }

  message("Python environment OK: ", py_config$python)
  invisible(TRUE)
}
