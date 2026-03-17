# Function to calculate the Affinity Matrix
affs <- function(x, na.action = na.omit) {
  x <- na.action(x)
  rown <- rownames(x)
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

# Estimate the optimal number of clusters (Rotation Cost Criterion)
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
    D = diag(degs)
    L = D - W
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
        sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
    }
    t2 <- sort(unlist(quality), index.return = TRUE)$ix
    K1 <- NUMC[t2[1]]
    K2 <- NUMC[t2[2]]
  }

  x <- list(first = K1, second = K2, value = unlist(quality))
  return(x)
}

# Internal utility: Discrete data processing for eigenvectors
.discretisationEigenVectorData <- function(eigenVector) {
  Y = matrix(0, nrow(eigenVector), ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector, 1, maxi)
  Y[cbind(1:nrow(eigenVector), j)] = 1
  return(Y)
}

# Internal utility: Spectral clustering discretization algorithm
.discretisation <- function(eigenVectors) {
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors, 1, normalize))

  n = nrow(eigenVectors)
  k = ncol(eigenVectors)

  R = matrix(0, k, k)
  R[, 1] = t(eigenVectors[round(n / 2), ])

  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }

  c = matrix(0, n, 1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[, j - 1], k, 1))
    i = mini(c)
    R[, j] = t(eigenVectors[i, ])
  }

  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]

    NcutValue = 2 * (n - sum(S))
    if (abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break

    lastObjectiveValue = NcutValue
    R = V %*% t(U)
  }

  return(list(discrete = eigenDiscrete, continuous = eigenVectors))
}

# Internal utility: Label propagation for prediction
.csPrediction <- function(W, Y0, method) {
  alpha = 0.9
  P = W / rowSums(W)
  if (method == 0) {
    Y = (1 - alpha) * solve(diag(dim(P)[1]) - alpha * P) %*% Y0
  } else {
    NLabel = which(rowSums(Y0) == 0)[1] - 1
    Y = Y0
    for (i in 1:1000) {
      Y = P %*% Y
      Y[1:NLabel, ] = Y0[1:NLabel, ]
    }
  }
  return(Y)
}

# Internal utility: Extract the top KK neighbors
.dominateset <- function(xx, KK = 20) {
  zero <- function(x) {
    s = sort(x, index.return = TRUE)
    x[s$ix[1:(length(x) - KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0, nrow(xx), ncol(xx))
  for (i in 1:nrow(xx)) {
    A[i, ] = zero(xx[i, ])
  }
  return(normalize(A))
}

# Calculate mutual information between vectors x and y
.mutualInformation <- function(x, y) {
  classx <- unique(x)
  classy <- unique(y)
  nx <- length(x)
  ncx <- length(classx)
  ncy <- length(classy)

  probxy <- matrix(NA, ncx, ncy)
  for (i in 1:ncx) {
    for (j in 1:ncy) {
      probxy[i, j] <- sum((x == classx[i]) & (y == classy[j])) / nx
    }
  }

  probx <- matrix(rowSums(probxy), ncx, ncy)
  proby <- matrix(colSums(probxy), ncx, ncy, byrow = TRUE)
  result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm = TRUE)
  return(result)
}

# Calculate entropy of vector x
.entropy <- function(x) {
  class <- unique(x)
  nx <- length(x)
  nc <- length(class)

  prob <- rep.int(NA, nc)
  for (i in 1:nc) {
    prob[i] <- sum(x == class[i]) / nx
  }

  result <- -sum(prob * log(prob, 2))
  return(result)
}

# R equivalent of MATLAB's repmat
.repmat = function(X, m, n) {
  if (is.null(dim(X))) {
    mx = length(X)
    nx = 1
  } else {
    mx = dim(X)[1]
    nx = dim(X)[2]
  }
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}

# Calculate proportions
perc <- function(x) {
  sapply(1:length(x), function(i) {
    x[i] / sum(x)
  })
}

#' Plot Kaplan-Meier Survival Curves
#'
#' @description Generates Kaplan-Meier survival curves based on clustering results.
#'
#' @param clinical Data frame containing survival data.
#' @param labels Vector of cluster labels.
#' @param limit Time limit for follow-up.
#' @param annot Annotation labels.
#' @param color Custom color palette.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Plot title.
#' @param legend.pos Legend position.
#' @param palette Color palette.
#' @param risk.table Whether to show risk table.
#' @param risk.table.ratio Ratio of risk table height.
#' @param anno.pos Position of p-value annotations.
#' @param anno.x.shift X-axis shift for p-value annotations.
#'
#' @import ggplot2
#' @import survival
#' @importFrom survminer ggsurvplot
#' @export
plot_KM <- function(clinical, labels, limit = NULL, annot = NULL, color = NULL,
                    xlab = "Follow up", ylab = "Survival Probability",
                    title = NULL, legend.pos = "top", palette = "jama_classic",
                    risk.table = T, risk.table.ratio = 0.4, anno.pos = "bottom",
                    anno.x.shift = 0.5) {

  # Main KM Plotting Function
  plot_KM <- function(clinical, labels, limit = NULL, annot = NULL, color = NULL,
                      xlab = "Follow up", ylab = "Survival Probability",
                      title = NULL, legend.pos = "top", palette = "jama_classic",
                      risk.table = T, risk.table.ratio = 0.4, anno.pos = "bottom",
                      anno.x.shift = 0.5) {

    # 1. Data Pre-processing
    time <- clinical[, 1]
    event <- clinical[, 2] == 1
    if (!is.null(limit)) {
      event[time > limit] <- FALSE
      time[time > limit] <- limit
    }
    df <- data.frame(futime = time, fustat = event, group = labels)

    # 2. Survival Analysis Calculation
    surv <- survival::survfit(survival::Surv(futime, fustat) ~ group, data = df)
    survstats <- survival::survdiff(survival::Surv(futime, fustat) ~ group, data = df)
    survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) - 1)

    # 3. Color Handling
    if (!is.null(color)) {
      if (!is.null(names(color))) {
        labels <- factor(labels, levels = names(color))
      }
    } else {
      color <- get_color(palette, n = length(unique(labels)))
    }

    # 4. Label Handling
    if (inherits(labels, "factor")) {
      legend.labs <- na.omit(levels(droplevels(labels[!(is.na(time) | is.na(event))])))
    } else if (is.logical(labels)) {
      labels <- factor(labels, levels = c(FALSE, TRUE))
      legend.labs <- na.omit(levels(droplevels(labels)))
    } else {
      legend.labs <- na.omit(unique(labels))
      labels <- factor(labels, levels = legend.labs)
    }

    # 5. Scientific Notation Helper
    fancy_scientific <- function(l, dig = 3) {
      l <- format(l, digits = dig, scientific = TRUE)
      l <- gsub("^(.*)e", "'\\1'e", l)
      l <- gsub("e", "%*%10^", l)
      parse(text = l)
    }

    # 6. Core Plotting
    p <- survminer::ggsurvplot(surv, data = df, xlab = xlab,
                               ylab = ylab, palette = color, legend = legend.pos,
                               legend.labs = legend.labs,
                               risk.table = risk.table,
                               risk.table.title = ggplot2::element_blank(),
                               risk.table.y.text = FALSE,
                               ggtheme = cowplot::theme_cowplot())

    p$plot <- p$plot + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.title = ggplot2::element_blank())

    # 7. Annotation Preparation
    anno.text <- ifelse(survstats$p.value == 0, "italic(P)<1%*%10^{-22}",
                        paste0("italic(P)==", fancy_scientific(survstats$p.value, 3)))
    anno.y.shift <- 0

    if (length(legend.labs) == 2) {
      hr <- survcomp::hazard.ratio(labels[!(is.na(time) | is.na(event))],
                                   time[!(is.na(time) | is.na(event))],
                                   event[!(is.na(time) | is.na(event))])
      anno.text <- c(anno.text, sprintf("HR == %3.2f~(%3.2f - %3.2f)",
                                        hr$hazard.ratio, hr$lower, hr$upper))
      anno.y.shift <- c(anno.y.shift + 0.15, 0)
    }

    if (!is.null(annot)) {
      anno.text <- c(anno.text, annot)
      anno.y.shift <- c(anno.y.shift + 0.15, 0)
    }

    # 8. Add Annotations to Plot
    if (anno.pos == "bottom") {
      p$plot <- p$plot + ggplot2::annotate("text", x = 0, y = anno.y.shift,
                                           label = anno.text, hjust = 0, vjust = 0, parse = TRUE)
    } else {
      p$plot <- p$plot + ggplot2::annotate("text", x = anno.x.shift * max(time, na.rm = TRUE),
                                           y = 0.85 + anno.y.shift, label = anno.text,
                                           hjust = 0, vjust = 2, parse = TRUE)
    }

    # 9. Final Composition with Risk Table
    if (risk.table) {
      p$table <- p$table + ggplot2::theme(axis.title.y = ggplot2::element_blank())
      pp <- cowplot::plot_grid(plotlist = list(p$plot + ggplot2::theme(axis.title.x = ggplot2::element_blank()),
                                               p$table + ggplot2::labs(x = xlab)),
                               labels = "", ncol = 1, align = "v",
                               rel_heights = c(1, risk.table.ratio))
      return(pp)
    } else {
      return(p$plot)
    }
  }

  # Color palette helper
  get_color <- function(palette, n = 12) {
    if(length(palette) > 1) return(palette)

    switch(tolower(palette), nature = {
      (ggsci::pal_npg("nrc"))(n)
    }, jco = {
      (ggsci::pal_jco("default"))(n)
    }, lancet = {
      (ggsci::pal_lancet("lanonc"))(n)
    }, jama = {
      ggsci::pal_jama()(n)
    }, jama_classic = {
      head(c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282","#FE850F","#E6580E","#FC491C","#FC491C","#E6190E","#FE0F5A","#E9FF63","#E6A00E"), n)
    },
    RColorBrewer::brewer.pal(n, "Set1")
    )
  }

  # Generate time-specific binary event matrix
  generate_time_event <- function(clinical, limits, labels = NULL) {
    time <- clinical[, 1]
    event <- clinical[, 2] == 1
    df <- sapply(limits, function(limit) {
      res <- event
      res[time > limit] <- F
      res
    })
    colnames(df) <- labels
    df
  }

  # Default theme settings
  .onLoad <- function(libname, pkgname) {
    ggplot2::theme_set(cowplot::theme_cowplot())
  }

  # Self-Diffusion algorithm for Affinity enhancement
  self.diffusion <- function(A, K) {
    diag(A) = 0
    sign_A = A
    sign_A[which(A > 0, arr.ind = TRUE)] = 1
    sign_A[which(A < 0, arr.ind = TRUE)] = -1

    P = dominate.set(abs(A), min(K, nrow(A) - 1)) * sign_A
    DD = apply(abs(P), MARGIN = 1, FUN = sum)
    diag(P) = DD + 1
    P = transition.fields(P)

    eigen_P = eigen(P)
    U = eigen_P$vectors
    D = eigen_P$values
    d = Re(D + .Machine$double.eps)

    alpha = 0.8
    beta = 2
    d = ((1 - alpha) * d) / (1 - alpha * d^beta)

    D = array(0, c(length(Re(d)), length(Re(d))))
    diag(D) = Re(d)

    W = U %*% D %*% t(U)
    diagonal_matrix = array(0, c(nrow(W), ncol(W)))
    diag(diagonal_matrix) = 1
    W = (W * (1 - diagonal_matrix)) / apply(array(0, c(nrow(W), ncol(W))), MARGIN = 2, FUN = function(x) {x = (1 - diag(W))})
    diag(D) = diag(D)[length(diag(D)):1]
    W = diag(DD) %*% W
    W = (W + t(W)) / 2
    W[which(W < 0, arr.ind = TRUE)] = 0

    return(W)
  }

  # Internal utility: Compute dominating set for K-nearest neighbors
  "dominate.set" <- function(aff.matrix, NR.OF.KNN) {
    PNN.matrix = array(0, c(nrow(aff.matrix), ncol(aff.matrix)))
    res.sort = apply(t(aff.matrix), MARGIN = 2, FUN = function(x) {return(sort(x, decreasing = TRUE, index.return = TRUE))})
    sorted.aff.matrix = t(apply(as.matrix(1:length(res.sort)), MARGIN = 1, function(x) { return(res.sort[[x]]$x) }))
    sorted.indices = t(apply(as.matrix(1:length(res.sort)), MARGIN = 1, function(x) { return(res.sort[[x]]$ix) }))

    res = sorted.aff.matrix[, 1:NR.OF.KNN]
    inds = array(0, c(nrow(aff.matrix), NR.OF.KNN))
    inds = apply(inds, MARGIN = 2, FUN = function(x) {x = 1:nrow(aff.matrix)})
    loc = sorted.indices[, 1:NR.OF.KNN]

    PNN.matrix[(as.vector(loc) - 1) * nrow(aff.matrix) + as.vector(inds)] = as.vector(res)
    PNN.matrix = (PNN.matrix + t(PNN.matrix)) / 2
    return(PNN.matrix)
  }

  # Internal utility: Compute transition fields for a matrix
  "transition.fields" <- function(W) {
    zero.index = which(apply(W, MARGIN = 1, FUN = sum) == 0)
    W = dn(W, 'ave')
    w = sqrt(apply(abs(W), MARGIN = 2, FUN = sum) + .Machine$double.eps)
    W = W / t(apply(array(0, c(nrow(W), ncol(W))), MARGIN = 2, FUN = function(x) {x = w}))
    W = W %*% t(W)
    W[zero.index, ] = 0
    W[, zero.index] = 0
    return(W)
  }

  # Normalization for symmetric kernels
  "dn" = function(w, type) {
    D = apply(w, MARGIN = 2, FUN = sum)
    if(type == "ave") {
      D = 1 / D
      D_temp = matrix(0, nrow = length(D), ncol = length(D))
      D_temp[cbind(1:length(D), 1:length(D))] = D
      D = D_temp
      wn = D %*% w
    } else if(type == "gph") {
      D = 1 / sqrt(D)
      D_temp = matrix(0, nrow = length(D), ncol = length(D))
      D_temp[cbind(1:length(D), 1:length(D))] = D
      D = D_temp
      wn = D %*% (w %*% D)
    } else {
      stop("Invalid type!")
    }
    return(wn)
  }

  # Core Spectral Clustering function
  spec.clu <- function(affinity, K, type = 3) {
    d = rowSums(affinity)
    d[d == 0] = .Machine$double.eps
    D = diag(d)
    L = D - affinity
    if (type == 1) {
      NL = L
    } else if (type == 2) {
      Di = diag(1 / d)
      NL = Di %*% L
    } else if(type == 3) {
      Di = diag(1 / sqrt(d))
      NL = Di %*% L %*% Di
    }
    eig = eigen(NL)
    res = sort(abs(eig$values), index.return = TRUE)
    U = eig$vectors[, res$ix[1:K]]
    normalize <- function(x) x / sqrt(sum(x^2))
    if (type == 3) {
      U = t(apply(U, 1, normalize))
    }
    eigDiscrete = .discretisation(U)
    eigDiscrete = eigDiscrete$discrete
    labels = apply(eigDiscrete, 1, which.max)
    return(labels)
  }
