#' @title Data Normalization Preprocessing (Internal)
#' @description Performs Z-score normalization on the input matrix and handles missing values.
#' @param data A numeric matrix where samples are rows and features are columns.
#' @return A standardized numeric matrix with NAs replaced by 0.
#' @keywords internal
preprocess_mydat <- function(data) {
  # Apply Z-score normalization column-wise: (x - mean) / sd
  data_scaled <- apply(data, 2, function(x) {
    s_dev <- sd(x, na.rm = TRUE)
    # Handle zero-variance features to avoid NaN
    if(is.na(s_dev) || s_dev == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s_dev
  })
  # Replace resulting NA values (e.g., from missing data) with 0
  data_scaled[is.na(data_scaled)] <- 0
  return(data_scaled)
}

#' @title OSSP Multi-omics Integrative Clustering Analysis
#' @description Executes the core pipeline including multi-omics data integration,
#' self-diffusion-based feature extraction, and K-means clustering.
#' @param data_list A list of numeric matrices (e.g., list(GE, ME, MI)).
#' @param K_clusters The number of clusters to be formed. Default is 7.
#' @param type The type of Laplacian matrix; 1: Unnormalized, 2: Random Walk, 3: Symmetric (Default).
#' @return A list containing:
#' \itemize{
#'   \item \code{labels}: A vector of integers indicating the cluster assignment for each sample.
#'   \item \code{reduced_data}: The integrated feature matrix after dimensionality reduction.
#' }
#' @importFrom stats kmeans
#' @export
run_ossp_analysis <- function(data_list, K_clusters = 7, type = 3) {

  # 1. Pre-processing
  message("Step 1: Pre-processing multi-omics data...")
  d <- lapply(data_list, preprocess_mydat)

  # 2. Eigenvector extraction
  message("Step 2: Extracting spectral eigenvectors...")
  ul <- lapply(seq_along(d), function(i) {
    # Calculate affinity and apply self-diffusion (assumes dependency functions are available)
    x <- affs(d[[i]])
    affinity <- self.diffusion(x, 4)

    # Construct Laplacian Matrix
    d_vec <- rowSums(affinity)
    # Use machine epsilon to avoid division by zero
    d_vec[d_vec == 0] <- .Machine$double.eps

    if (type == 1) {
      NL <- diag(d_vec) - affinity
    } else if (type == 2) {
      NL <- diag(length(d_vec)) - diag(1 / d_vec) %*% affinity
    } else {
      # Symmetric Normalized Laplacian
      Di <- diag(1 / sqrt(d_vec))
      NL <- diag(length(d_vec)) - Di %*% affinity %*% Di
    }

    eig <- eigen(NL)
    res <- sort(abs(eig$values), index.return = TRUE)

    # Automated feature dimension determination (Eigengap heuristic)
    e <- res$x[1:15]
    sa <- sapply(1:14, function(j) abs((e[j+1] - e[j])))
    K_auto <- which(sa == max(sa)) + 3

    # Extract and normalize eigenvectors
    U <- eig$vectors[, res$ix[1:K_auto]]
    if (type == 3) {
      U <- t(apply(U, 1, function(x) x / sqrt(sum(x^2))))
    }
    return(U)
  })

  # 3. Concatenation and Clustering
  message("Step 3: Integrating features and performing K-means clustering...")
  rd <- do.call(cbind, ul)
  set.seed(11111)
  labx <- kmeans(rd, centers = K_clusters)

  return(list(labels = labx$cluster, reduced_data = rd))
}

#' @title Evaluate Clustering Quality (Silhouette Analysis)
#' @description Computes the fused similarity matrix and generates a Silhouette plot
#' to assess the consistency within clusters.
#' @param reduced_data The fused feature matrix returned by \code{run_ossp_analysis}.
#' @param labels A vector of cluster labels.
#' @return A silhouette object containing the calculated widths.
#' @export
eval_ossp_clustering <- function(reduced_data, labels) {
  message("Step 4: Evaluating clustering quality...")
  # Compute similarity matrix for the integrated features
  x_sim <- affs(reduced_data)

  # Calculate Silhouette widths (requires silhouette_SimilarityMatrix function)
  sil <- silhouette_SimilarityMatrix(labels, x_sim)

  # Generate Plot
  plot(sil, main = "OSSP Clustering Silhouette Plot")
  return(sil)
}

#' @title Plot Kaplan-Meier Survival Curves
#' @description Generates Kaplan-Meier survival plots based on OSSP clustering results.
#' @param survival_data A data frame containing survival metadata with 'Survival' and 'Death' columns.
#' @param labels A vector of cluster labels.
#' @importFrom survival Surv
#' @return A survival plot object.
#' @export
plot_ossp_km <- function(survival_data, labels) {
  message("Step 5: Plotting Kaplan-Meier survival curves...")
  time <- as.numeric(survival_data$Survival)
  event <- as.numeric(survival_data$Death)
  clins <- survival::Surv(time, event)

  # Call KM plotting function
  plot_KM(
    clins,
    as.integer(labels),
    palette = "lanonc",
    xlab = "Follow up (Months)",
    ylab = "Overall Survival Probability",
    risk.table = TRUE
  )
}
