# 函数 1：数据标准化
#' 数据标准化与缺失值处理
#'
#' @description 对输入的矩阵进行 Z-score 标准化，并将产生的 NA 值填充为 0。
#'
#' @param data 数值矩阵（样本为行，特征为列）。
#'
#' @return 返回经过标准化处理后的数值矩阵。
preprocess_mydat <- function(data) {
  data_scaled <- apply(data, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  data_scaled[is.na(data_scaled)] <- 0
  return(data_scaled)
}

# 函数 2：主分析程序
#' OSSP 多组学整合聚类分析
#'
#' @description 执行多组学数据的整合、特征提取及 K-means 聚类。
#'
#' @param mydatGE 基因表达谱矩阵（样本行，基因列）。
#' @param mydatME DNA 甲基化数据矩阵（样本行，探针列）。
#' @param mydatMI miRNA 表达数据矩阵（样本行，miRNA 列）。
#' @param K_clusters 预设聚类数，默认为 7。
#'
#' @return 包含聚类标签 (labels) 和融合特征矩阵 (reduced_data) 的列表。
#' @export
run_ossp_analysis <- function(mydatGE, mydatME, mydatMI, K_clusters = 7) {
  # 预处理
  d <- list(
    preprocess_mydat(mydatGE),
    preprocess_mydat(mydatME),
    preprocess_mydat(mydatMI)
  )

  # 特征向量提取逻辑
  ul <- lapply(1:3, function(i) {
    x <- affs(d[[i]])
    affinity <- self.diffusion(x, 4)

    # 谱聚类核心逻辑
    d_vec <- rowSums(affinity)
    d_vec[d_vec == 0] <- .Machine$double.eps
    Di <- diag(1 / sqrt(d_vec))
    L <- diag(d_vec) - affinity
    NL <- Di %*% L %*% Di

    eig <- eigen(NL)
    res <- sort(abs(eig$values), index.return = TRUE)

    # 自动判定特征维度
    e <- res$x[1:15]
    sa <- sapply(1:14, function(i) abs((e[i+1] - e[i])))
    K_auto <- which(sa == max(sa)) + 3

    U <- eig$vectors[, res$ix[1:K_auto]]
    U <- t(apply(U, 1, function(x) x / sqrt(sum(x^2))))
    return(U)
  })

  # 合并与最终聚类
  rd <- do.call(cbind, ul)
  set.seed(11111)
  labx <- kmeans(rd, K_clusters)

  return(list(labels = labx$cluster, reduced_data = rd))
}

# 函数 3：结果可视化
#' 绘制聚类生存分析 KM 曲线
#'
#' @description 基于聚类结果绘制 Kaplan-Meier 生存曲线。
#'
#' @param survival_data 生存数据框，需包含 'Survival' 和 'Death' 列。
#' @param labels 聚类标签向量。
#'
#' @import ggplot2
#' @importFrom survival Surv
#' @return 绘制并显示生存曲线图。
#' @export
plot_ossp_km <- function(survival_data, labels) {
  time <- as.numeric(survival_data$Survival)
  event <- as.numeric(survival_data$Death)
  clins <- cbind(time, event) 
  plot_KM(
    clins,
    as.integer(labels),
    palette = "lanonc",
    xlab = "Follow up(Months)",
    ylab = "OS(pro.)",
    risk.table = TRUE
  )
}
