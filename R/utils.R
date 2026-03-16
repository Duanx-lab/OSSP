affs <- function(x, na.action = na.omit) {
  # 处理缺失值
  x <- na.action(x)
  rown <- rownames(x)
  x <- as.matrix(x)
  m <- nrow(x)

  # 初始化向量和矩阵
  s <- rep(0, m)
  dota <- rowSums(x * x) / 2
  dis <- crossprod(t(x))  # dis 是一个 m x m 的矩阵



  # 计算距离矩阵
  for (i in 1:m) {
    dis[i, ] <- 2 * (-dis[i, ] + dota + rep(dota[i], m))
  }

  # 防止负数出现
  dis[dis < 0] <- 0

  # 计算每行的中位数
  for (i in 1:m) {
    s[i] <- median(sort(sqrt(dis[i, ]))[1:5])
  }



  # 计算 Affinity 矩阵
  km <- exp(-dis / (s %*% t(s)))  # 确保维度正确，避免广播错误

  return(km)
}


# 估算图的最佳聚类数 (Rotation Cost 准则)
# 供内部调用，评估最佳 (first) 和次优 (second) 聚类数
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

# 内部辅助：特征向量离散化数据处理
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

# 内部辅助：谱聚类离散化核心算法
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



.csPrediction <- function(W,Y0,method){
  ###This function implements the label propagation to predict the label(subtype) for new patients.
  ### note method is an indicator of which semi-supervised method to use
  # method == 0 indicates to use the local and global consistency method
  # method >0 indicates to use label propagation method.

  alpha=0.9;
  P= W/rowSums(W)
  if(method==0){
    Y= (1-alpha)* solve( diag(dim(P)[1])- alpha*P)%*%Y0;
  } else {
    NLabel=which(rowSums(Y0)==0)[1]-1;
    Y=Y0;
    for (i in 1:1000){
      Y=P%*%Y;
      Y[1:NLabel,]=Y0[1:NLabel,];
    }
  }
  return(Y);
}

.discretisation <- function(eigenVectors) {

  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))

  n = nrow(eigenVectors)
  k = ncol(eigenVectors)

  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])

  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }

  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }

  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)

    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]

    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break

    lastObjectiveValue = NcutValue
    R = V %*% t(U)

  }

  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {

  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1

  return(Y)

}

.dominateset <- function(xx,KK=20) {
  ###This function outputs the top KK neighbors.

  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
  }
  return(normalize(A))
}

# Calculate the mutual information between vectors x and y.
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
  proby <- matrix(colSums(probxy), ncx, ncy, byrow=TRUE)
  result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm=TRUE)
  return(result)
}

# Calculate the entropy of vector x.
.entropy <- function(x) {
  class <- unique(x)
  nx <- length(x)
  nc <- length(class)

  prob <- rep.int(NA, nc)
  for (i in 1:nc) {
    prob[i] <- sum(x == class[i])/nx
  }

  result <- -sum(prob * log(prob, 2))
  return(result)
}

.repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  if (is.null(dim(X))) {
    mx = length(X)
    nx = 1
  } else {
    mx = dim(X)[1]
    nx = dim(X)[2]
  }
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}


perc<-function(x){
  sapply(1:length(x), function(i){
    x[i]/sum(x)
  })
}

#' @import ggplot2
#' @import survival
#' @importFrom survminer ggsurvplot
#' @importFrom cowplot theme_cowplot plot_grid
#' @importFrom ggsci pal_npg pal_jco pal_lancet pal_jama
#' @importFrom survcomp hazard.ratio
#' @importFrom stats pchisq na.omit
#' @importFrom utils head
#' @export
plot_KM <- function (clinical, labels, limit = NULL, annot = NULL, color = NULL,
                     xlab = "Follow up", ylab = "Survival Probability",
                     title = NULL, legend.pos = "top", palette = "jama_classic",
                     risk.table = T, risk.table.ratio = 0.4, anno.pos = "bottom",
                     anno.x.shift = 0.5)
{
  # ... 你原来的代码函数体 ...
}
{
  time <- clinical[, 1]
  event <- clinical[, 2] == 1
  if (!is.null(limit)) {
    event[time > limit] <- F
    time[time > limit] <- limit
  }
  df <- data.frame(futime = time, fustat = event, group = labels)
  surv <- survival::survfit(survival::Surv(futime, fustat) ~
                              group, data = df)
  survstats <- survival::survdiff(survival::Surv(futime, fustat) ~
                                    group, data = df)
  survstats$p.value <- 1 - pchisq(survstats$chisq, length(survstats$n) -
                                    1)
  if (!is.null(color)) {
    if (!is.null(names(color))) {
      labels <- factor(labels, levels = names(color))
    }
  }
  else {
    color <- get_color(palette, n = length(unique(labels)))
  }
  if (class(labels) == "factor") {
    legend.labs <- na.omit(levels(droplevels(labels[!(is.na(time) |
                                                        is.na(event))])))
  }
  else if (class(labels) == "logical") {
    labels <- factor(labels, levels = c(F, T))
    legend.labs <- na.omit(levels(droplevels(labels)))
  }
  else {
    legend.labs <- na.omit(unique(labels))
    labels <- factor(labels, levels = legend.labs)
  }
  fancy_scientific <- function(l, dig = 3) {
    l <- format(l, digits = dig, scientific = TRUE)
    l <- gsub("^(.*)e", "'\\1'e", l)
    l <- gsub("e", "%*%10^", l)
    parse(text = l)
  }
  p <- survminer::ggsurvplot(surv, data = df, xlab = xlab,
                             ylab = ylab, palette = color, legend = legend.pos, legend.labs = legend.labs,
                             risk.table = risk.table, risk.table.title = element_blank(),
                             risk.table.y.text = FALSE, ggtheme =cowplot::theme_cowplot())
  p$plot <- p$plot + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),

          legend.title = element_blank())
  anno.text <- ifelse(survstats$p.value == 0, "italic(P)<1%*%10^{-22}",
                      paste0("italic(P)==", fancy_scientific(survstats$p.value,
                                                             3)))
  anno.y.shift <- 0
  if (length(legend.labs) == 2) {
    hr <- survcomp::hazard.ratio(labels[!(is.na(time) | is.na(event))],
                                 time[!(is.na(time) | is.na(event))], event[!(is.na(time) |
                                                                                is.na(event))])
    anno.text <- c(anno.text, sprintf("HR == %3.2f~(%3.2f - %3.2f)",
                                      hr$hazard.ratio, hr$lower, hr$upper))
    anno.y.shift <- c(anno.y.shift + 0.15, 0)
  }
  if (!is.null(annot)) {
    anno.text <- c(anno.text, annot)
    anno.y.shift <- c(anno.y.shift + 0.15, 0)
  }
  if (anno.pos == "bottom") {
    p$plot <- p$plot + annotate("text", x = 0,
                                y = anno.y.shift, label = anno.text, hjust = 0, vjust = 0,
                                parse = TRUE)
  }
  else {
    p$plot <- p$plot + annotate("text", x = anno.x.shift *
                                  max(time, na.rm = T), y = 0.85 + anno.y.shift, label = anno.text,
                                hjust = 0, vjust = 2, parse = TRUE)
  }
  if (risk.table) {
    p$table <- p$table  +  theme(

      axis.title.y = element_blank())
    pp <- cowplot::plot_grid(plotlist = list(p$plot + theme(axis.title.x = element_blank()),
                                             p$table + labs(x = xlab)), labels = "", ncol = 1,
                             align = "v", rel_heights = c(1, risk.table.ratio))
    return(pp)
  }
  else return(p$plot)
}



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

#' @export
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


.onLoad <- function(libname, pkgname) {
  ggplot2::theme_set(cowplot::theme_cowplot())
}



get_color <- function(palette, n = 6) {
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
    head(c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"), n)
  },
  RColorBrewer::brewer.pal(n, "Set1")
  )
}




self.diffusion <- function(A,K) {

  # set the values of the diagonal of A to 0
  diag(A) = 0

  # compute the sign matrix of A
  sign_A = A
  sign_A[which(A>0,arr.ind=TRUE)] = 1
  sign_A[which(A<0,arr.ind=TRUE)] = -1

  # compute the dominate set for A and K
  P = dominate.set(abs(A),min(K,nrow(A)-1)) * sign_A

  # sum the absolute value of each row of P
  DD = apply(abs(P),MARGIN=1,FUN=sum)

  # set DD+1 to the diagonal of P
  diag(P) = DD + 1

  # compute the transition field of P
  P = transition.fields(P)

  # compute the eigenvalues and eigenvectors of P
  eigen_P = eigen(P)
  U = eigen_P$vectors
  D = eigen_P$values

  # set to d the real part of the diagonal of D
  d = Re(D + .Machine$double.eps)

  #perform the diffusion
  alpha = 0.8
  beta = 2
  d = ((1-alpha)*d)/(1-alpha*d^beta)

  # set to D the real part of the diagonal of d
  D = array(0,c(length(Re(d)),length(Re(d))))
  diag(D) = Re(d)

  # finally compute W
  W = U %*% D %*% t(U)
  diagonal_matrix = array(0,c(nrow(W),ncol(W)))
  diag(diagonal_matrix) = 1
  W = (W * (1-diagonal_matrix)) / apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=(1-diag(W))})
  diag(D) = diag(D)[length(diag(D)):1]
  W = diag(DD) %*% W
  W = (W + t(W)) / 2

  W[which(W<0,arr.ind=TRUE)] = 0

  return(W)

}

# compute the dominate set for the matrix aff.matrix and NR.OF.KNN
"dominate.set" <- function( aff.matrix, NR.OF.KNN ) {

  # create the structure to save the results
  PNN.matrix = array(0,c(nrow(aff.matrix),ncol(aff.matrix)))

  # sort each row of aff.matrix in descending order and saves the sorted
  # array and a collection of vectors with the original indices
  res.sort = apply(t(aff.matrix),MARGIN=2,FUN=function(x) {return(sort(x, decreasing = TRUE, index.return = TRUE))})
  sorted.aff.matrix = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$x) }))
  sorted.indices = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$ix) }))

  # get the first NR.OF.KNN columns of the sorted array
  res = sorted.aff.matrix[,1:NR.OF.KNN]

  # create a matrix of NR.OF.KNN columns by binding vectors of
  # integers from 1 to the number of rows/columns of aff.matrix
  inds = array(0,c(nrow(aff.matrix),NR.OF.KNN))
  inds = apply(inds,MARGIN=2,FUN=function(x) {x=1:nrow(aff.matrix)})

  # get the first NR.OF.KNN columns of the indices of aff.matrix
  loc = sorted.indices[,1:NR.OF.KNN]

  # assign to PNN.matrix the sorted indices
  PNN.matrix[(as.vector(loc)-1)*nrow(aff.matrix)+as.vector(inds)] = as.vector(res)

  # compute the final results and return them
  PNN.matrix = (PNN.matrix + t(PNN.matrix))/2

  return(PNN.matrix)

}

# compute the transition field of the given matrix
"transition.fields" <- function( W ) {

  # get any index of columns with all 0s
  zero.index = which(apply(W,MARGIN=1,FUN=sum)==0)

  # compute the transition fields
  W = dn(W,'ave')

  w = sqrt(apply(abs(W),MARGIN=2,FUN=sum)+.Machine$double.eps)
  W = W / t(apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=w}))
  W = W %*% t(W)

  # set to 0 the elements of zero.index
  W[zero.index,] = 0
  W[,zero.index] = 0

  return(W)

}

# normalizes a symmetric kernel
"dn" = function( w, type ) {

  # compute the sum of any column
  D = apply(w,MARGIN=2,FUN=sum)

  # type "ave" returns D^-1*W
  if(type=="ave") {
    D = 1 / D
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% w
  }
  # type "gph" returns D^-1/2*W*D^-1/2
  else if(type=="gph") {
    D = 1 / sqrt(D)
    D_temp = matrix(0,nrow=length(D),ncol=length(D))
    D_temp[cbind(1:length(D),1:length(D))] = D
    D = D_temp
    wn = D %*% (w %*% D)
  }
  else {
    stop("Invalid type!")
  }

  return(wn)

}


spec.clu<- function(affinity, K, type=3) {


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
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)



  return(labels)
}




.discretisation <- function(eigenVectors) {

  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))

  n = nrow(eigenVectors)
  k = ncol(eigenVectors)

  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])

  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }

  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }

  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)

    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]

    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps)
      break

    lastObjectiveValue = NcutValue
    R = V %*% t(U)

  }

  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}
