#' Run a complete spatial metabolomics clustering workflow
#'
#' This high-level function performs a complete workflow for spatial metabolomics
#' clustering:
#' \enumerate{
#'   \item Extract spectra matrix and pixel metadata from a Cardinal MSI object
#'   \item Remove constant features
#'   \item Apply min-max normalization
#'   \item Run Python UMAP through \pkg{reticulate}
#'   \item Perform k-means clustering on UMAP coordinates
#'   \item Build a clustering result table
#'   \item Write cluster labels back into \code{pixelData(msi_obj)}
#' }
#'
#' This workflow requires a Python environment configured through
#' \pkg{reticulate} when UMAP is run via the Python backend.
#'
#' @param msi_obj A Cardinal MSI object.
#' @param python_path Optional character string specifying the Python executable.
#' If \code{NULL}, the current \pkg{reticulate} Python configuration will be used.
#' @param centers Integer; number of k-means clusters. Default is \code{10}.
#' @param metric Character; UMAP distance metric. Default is \code{"cosine"}.
#' @param n_neighbors Integer; number of neighbors for UMAP. Default is \code{10L}.
#' @param min_dist Numeric; UMAP \code{min_dist} parameter. Default is \code{0.05}.
#' @param n_components Integer; number of UMAP dimensions. Default is \code{2L}.
#' @param umap_seed Integer; random seed for UMAP. Default is \code{2025L}.
#' @param kmeans_seed Integer; random seed for k-means. Default is \code{2024}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{spectra_filtered}{The spectra matrix after removing constant features.}
#'   \item{spectra_scaled}{The min-max normalized spectra matrix.}
#'   \item{umap_df}{A data.frame containing UMAP coordinates.}
#'   \item{kmeans_result}{A list returned by \code{run_kmeans_cluster()},
#'   containing the raw \code{kmeans} result and cluster assignments.}
#'   \item{cluster_df}{A data.frame containing embedding, cluster labels, and metadata.}
#'   \item{msi_obj}{The updated Cardinal MSI object with cluster labels in pixelData.}
#' }
#'
#' @examples
#' \dontrun{
#' library(Cardinal)
#' msi_obj <- readImzML("example.imzML")
#'
#' res <- spatial_kmeans_workflow(
#'   msi_obj = msi_obj,
#'   python_path = "/path/to/python",
#'   centers = 10
#' )
#'
#' head(res$umap_df)
#' table(res$cluster_df$kmeans_cluster)
#' }
#'
#' @export
spatial_kmeans_workflow <- function(
    msi_obj,
    python_path = NULL,
    centers = 10,
    metric = "cosine",
    n_neighbors = 10L,
    min_dist = 0.05,
    n_components = 2L,
    umap_seed = 2025L,
    kmeans_seed = 2024
) {
  extracted <- extract_spectra_matrix(msi_obj)

  spectra <- extracted$spectra
  pixel_info <- extracted$pixel_info

  spectra <- remove_constant_features(spectra)
  spectra_scaled <- minmax_normalize(spectra)

  umap_df <- run_umap_py(
    x = spectra_scaled,
    python_path = python_path,
    metric = metric,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    n_components = n_components,
    random_state = umap_seed
  )

  km_res <- run_kmeans_cluster(
    embedding = umap_df,
    centers = centers,
    seed = kmeans_seed
  )

  cluster_df <- build_cluster_dataframe(
    embedding = umap_df,
    cluster = km_res$cluster,
    pixel_info = pixel_info
  )

  msi_obj_updated <- attach_cluster_to_pixeldata(
    msi_obj = msi_obj,
    cluster_df = cluster_df
  )

  list(
    spectra_filtered = spectra,
    spectra_scaled = spectra_scaled,
    umap_df = umap_df,
    kmeans_result = km_res,
    cluster_df = cluster_df,
    msi_obj = msi_obj_updated
  )

}

