#' @title Compute Affinity Matrix
#' @description Calculates the affinity matrix using Euclidean distances with local scaling.
#' @param x A numeric matrix (samples as rows, features as columns).
#' @param na.action Action to handle missing values. Defaults to \code{na.omit}.
#' @return A symmetric affinity matrix.
#' @export
affs <- function(x, na.action = na.omit) {
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)

  s <- rep(0, m)
  dota <- rowSums(x * x) / 2
  dis <- crossprod(t(x))

  for (i in 1:m) {
    dis[i, ] <- 2 * (-dis[i, ] + dota + rep(dota[i], m))
  }

  dis[dis < 0] <- 0
  for (i in 1:m) {
    s[i] <- median(sort(sqrt(dis[i, ]))[1:5])
  }

  km <- exp(-dis / (s %*% t(s)))
  return(km)
}

#' @title Estimate Optimal Number of Clusters
#' @description Estimates the optimal number of clusters for a given graph using rotation cost.
#' @param W The affinity graph matrix.
#' @param NUMC A vector of possible numbers of clusters. Defaults to 2:8.
#' @return A list containing the first best (\code{first}), second best (\code{second}), and all quality values.
#' @export
es.num.graph <- function(W, NUMC = 2:8) {
  if (min(NUMC) == 1) {
    warning('Note that we always assume there are more than one cluster.')
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W)) / 2
  diag(W) = 0

  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    L = diag(degs) - W
    Di = diag(1 / sqrt(degs))
    L = Di %*% L %*% Di
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = TRUE)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]

    quality = list()
    for (c_index in 1:length(NUMC)) {
      ck = NUMC[c_index]
      UU = eigs$vectors[, 1:ck]
      EigenvectorsDiscrete <- .discretisation(UU)[[1]]
      EigenVectors = EigenvectorsDiscrete^2
      temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
      temp1 <- t(apply(temp1, 1, sort, TRUE))

      quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) *
        sum(sum(diag(1 / (temp1[, 1] + .Machine$double.eps)) %*% temp1[, 1:max(2, ck - 1)]))
    }
    t2 <- sort(unlist(quality), index.return = TRUE)$ix
    K1 <- NUMC[t2[1]]
    K2 <- NUMC[t2[2]]
  }

  res <- list(first = K1, second = K2, value = unlist(quality))
  return(res)
}

#' @title Self-Diffusion on Affinity Graph
#' @description Implements self-diffusion logic to enhance local scaling affinity.
#' @param A Affinity graph matrix.
#' @param K Number of nearest neighbors for the dominate set.
#' @return A diffused symmetric matrix.
#' @export
self.diffusion <- function(A, K) {
  diag(A) = 0
  sign_A = A
  sign_A[A > 0] = 1
  sign_A[A < 0] = -1

  P = .dominateset(abs(A), min(K, nrow(A) - 1)) * sign_A
  DD = apply(abs(P), 1, sum)
  diag(P) = DD + 1

  P = .transition.fields(P)
  eigen_P = eigen(P)
  U = eigen_P$vectors
  d = Re(eigen_P$values + .Machine$double.eps)

  alpha = 0.8
  beta = 2
  d = ((1 - alpha) * d) / (1 - alpha * d^beta)

  D_mat = diag(Re(d))
  W = U %*% D_mat %*% t(U)

  diag_W = diag(W)
  W = (W * (1 - diag(nrow(W)))) / (1 - diag_W)
  W = diag(DD) %*% W
  W = (W + t(W)) / 2
  W[W < 0] = 0

  return(W)
}

#' @title Spectral Clustering
#' @description Performs spectral clustering on an affinity matrix.
#' @param affinity The input affinity graph.
#' @param K Number of clusters.
#' @param type Type of Laplacian; 1: Unnormalized, 2: Random Walk, 3: Symmetric (Default).
#' @return A vector of cluster labels.
#' @export
spec.clu <- function(affinity, K, type = 3) {
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  L = diag(d) - affinity

  if (type == 2) {
    L = diag(1 / d) %*% L
  } else if (type == 3) {
    Di = diag(1 / sqrt(d))
    L = Di %*% L %*% Di
  }

  eig = eigen(L)
  res = sort(abs(eig$values), index.return = TRUE)
  U = eig$vectors[, res$ix[1:K]]

  if (type == 3) {
    U = t(apply(U, 1, function(x) x / sqrt(sum(x^2))))
  }

  labels = apply(.discretisation(U)$discrete, 1, which.max)
  return(labels)
}

#' @title Plot Kaplan-Meier Curves
#' @description Generates survival plots with risk tables.
#' @param clinical A \code{Surv} object from the \code{survival} package.
#' @param labels A vector of cluster/group labels.
#' @param palette Color palette name. Defaults to "jama_classic".
#' @param risk.table Logical; whether to display a risk table.
#' @param ... Additional arguments passed to ggsurvplot.
#' @return A ggplot object (cowplot grid).
#' @import ggplot2
#' @importFrom survival survfit survdiff
#' @importFrom survminer ggsurvplot
#' @importFrom cowplot plot_grid
#' @export
plot_KM <- function(clinical, labels, palette = "jama_classic", risk.table = TRUE, ...) {
  # 内部逻辑保持不变
  df <- data.frame(futime = clinical[, 1], fustat = clinical[, 2], group = labels)
  surv_fit <- survival::survfit(survival::Surv(futime, fustat) ~ group, data = df)

  color_vals <- .get_color(palette, length(unique(labels)))

  p <- survminer::ggsurvplot(surv_fit, data = df, palette = color_vals,
                             risk.table = risk.table, ggtheme = cowplot::theme_cowplot(), ...)

  if (risk.table) {
    return(cowplot::plot_grid(p$plot, p$table, ncol = 1, align = "v", rel_heights = c(1, 0.4)))
  } else {
    return(p$plot)
  }
}

#' @title Generate Time-Event Matrix
#' @description Creates binary event indicators at specified time limits.
#' @param clinical A survival matrix or data frame.
#' @param limits A vector of time limits.
#' @param labels Labels for the output columns.
#' @export
generate_time_event <- function(clinical, limits, labels = NULL) {
  time <- clinical[, 1]
  event <- clinical[, 2] == 1
  df <- sapply(limits, function(limit) {
    res <- event
    res[time > limit] <- FALSE
    res
  })
  if (!is.null(labels)) colnames(df) <- labels
  return(df)
}

# --- Internal Helper Functions (Non-Exported) ---

.discretisation <- function(eigenVectors) {
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors, 1, normalize))
  n = nrow(eigenVectors); k = ncol(eigenVectors)
  R = matrix(0, k, k); R[, 1] = t(eigenVectors[round(n / 2), ])
  mini <- function(x) which(x == min(x))[1]
  c = matrix(0, n, 1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[, j - 1], k, 1))
    R[, j] = t(eigenVectors[mini(c), ])
  }
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    NcutValue = 2 * (n - sum(svde$d))
    if (abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) break
    lastObjectiveValue = NcutValue
    R = svde$v %*% t(svde$u)
  }
  return(list(discrete = eigenDiscrete, continuous = eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  Y = matrix(0, nrow(eigenVector), ncol(eigenVector))
  j = apply(eigenVector, 1, function(x) which(x == max(x))[1])
  Y[cbind(1:nrow(eigenVector), j)] = 1
  return(Y)
}

.dominateset <- function(xx, KK = 20) {
  zero <- function(x) {
    s = sort(x, index.return = TRUE)
    x[s$ix[1:(length(x) - KK)]] = 0
    return(x)
  }
  A = matrix(0, nrow(xx), ncol(xx))
  for (i in 1:nrow(xx)) A[i, ] = zero(xx[i, ])
  return(A / rowSums(A))
}

.transition.fields <- function(W) {
  zero.index = which(rowSums(W) == 0)
  W = .dn(W, 'ave')
  w = sqrt(colSums(abs(W)) + .Machine$double.eps)
  W = W / t(matrix(w, nrow(W), ncol(W), byrow = TRUE))
  W = W %*% t(W)
  W[zero.index, ] = 0; W[, zero.index] = 0
  return(W)
}

.dn <- function(w, type) {
  D = colSums(w)
  if (type == "ave") {
    wn = diag(1 / D) %*% w
  } else if (type == "gph") {
    Di = diag(1 / sqrt(D))
    wn = Di %*% w %*% Di
  }
  return(wn)
}

.get_color <- function(palette, n = 6) {
  jama_classic <- c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282")
  if (palette == "jama_classic") return(head(jama_classic, n))
  # Fallback to Set1 or expand via switch as needed
  return(RColorBrewer::brewer.pal(min(n, 9), "Set1"))
}


#' @title Evaluate Clustering Quality
#' @description Calculates the silhouette width and generates a silhouette plot based on the affinity matrix to assess clustering performance.
#' @param reduced_data The fused feature matrix (typically the output from \code{run_ossp_analysis}).
#' @param labels The cluster labels assigned to samples.
#' @return An object of class \code{silhouette} containing the calculated widths.
#' @importFrom cluster silhouette
#' @export
eval_ossp_clustering <- function(reduced_data, labels) {
  # 1. Compute the affinity matrix for the reduced feature space
  message("Calculating affinity matrix for evaluation...")
  x_sim <- affs(reduced_data)

  # 2. Compute the Silhouette coefficient
  # Convert similarity to distance (Distance = 1 - Similarity)
  # Ensure labels are treated as an integer vector for the silhouette function
  dist_matrix <- as.dist(1 - x_sim)
  sil <- cluster::silhouette(as.integer(as.factor(labels)), dist_matrix)

  # 3. Print summary statistics and generate the plot
  # Displays the average silhouette width and cluster-specific metrics in the console
  print(summary(sil))

  # Generate the silhouette diagram using the defined color palette
  plot(sil,
       main = "OSSP Clustering Silhouette Plot",
       col = .get_color("jama_classic", length(unique(labels))))

  return(sil)
}
