#set page(margin: 2.2cm)
#set text(size: 11pt, lang: "zh")
#set par(justify: true)
#set heading(numbering: "1.")

= 分段回归与季节分解模块

本模块位于 `R/Component.R`，实现季节均值分解与分段线性回归，供水文气象
时间序列分析使用。模块包含四个函数：

+ `RemoveSeasonal()`：季节均值与去季节序列；
+ `LS_OneBreak()`：单断点分段回归（同斜率、不同截距）；
+ `LS_MultiBreaks()`：多断点分段回归（同斜率、不同截距）；
+ `lm_solve()`：基于 `.lm.fit()` 的最小二乘内部求解器。

== RemoveSeasonal

#table(
  columns: (auto, 1fr),
  [参数], [说明],
  [`y`], [观测值向量。],
  [`t_mmdd`], [各观测对应的月日编码（MMDD），决定季节分组。],
)

返回值：

+ `mu`：按月日升序排列的季节均值，长度等于 `t_mmdd` 的唯一值个数；
+ `anomaly`：去季节序列，长度与 `y` 相同。

季节均值与异常值的定义为

$ mu_j = 1 / n_j sum_(t_i = j) y_i, quad
  "anomaly"_i = y_i - mu_(t_i), $

其中 $n_j$ 为第 $j$ 个季节组的样本数。均值计算忽略缺失值（`na.rm = TRUE`）。

```r
RemoveSeasonal <- function(y, t_mmdd) {
  inds <- sort(unique(t_mmdd))
  mu <- numeric(length(inds))
  anomaly <- numeric(length(y))
  for (i in seq_along(inds)) {
    ind <- which(t_mmdd == inds[i])
    mu[i] <- mean(y[ind], na.rm = TRUE)
    anomaly[ind] <- y[ind] - mu[i]
  }
  list(mu = mu, anomaly = anomaly)
}
```

== LS_OneBreak

单断点分段回归：断点前后共享斜率 $beta_1$，截距在断点 $I_c$ 处发生跳变。
模型写作

$ y_t = beta_0 + beta_1 x_t + beta_2 bold(1)_{t > I_c} + epsilon_t, $

其中 $bold(1)_{t > I_c}$ 为指示变量：$t <= I_c$ 时取 0，$t > I_c$ 时取 1。
设计矩阵为 $X = [bold(1), x, bold(1)_{t > I_c}]$。

#table(
  columns: (auto, 1fr),
  [参数], [说明],
  [`Y`], [响应值向量。],
  [`X`], [自变量向量。],
  [`Ic`], [变点位置：$[1, I_c]$ 与 $[I_c + 1, n]$ 分别对应前后两段。],
  [`...`], [透传给 `lm_solve()`。],
)

返回值：

+ `coef`：回归系数 $(beta_0, beta_1, beta_2)$；
+ `fitted`：拟合值；
+ `resi`：残差；
+ `SSE`：残差平方和。

== LS_MultiBreaks

多断点分段回归：$k$ 个断点将序列划分为 $k + 1$ 段，各段斜率相同、截距
不同。设 $I p s = (I_1, dots, I_k, n)$ 为各段末尾位置（末位为序列长度），
模型为

$ y_t = beta_0 + beta_1 (x_t - bar x) + sum_(j = 1)^k
  beta_(j + 1) bold(1)_{I_j < t <= I_(j + 1)} + epsilon_t . $

首段 $(1, I_1]$ 作为基准截距，其余 $k$ 段各引入一个指示变量。

#table(
  columns: (auto, 1fr),
  [参数], [说明],
  [`Y`], [响应值向量。],
  [`X`], [自变量向量。],
  [`Ips`], [各段末尾位置，末位为序列长度。],
  [`...`], [透传给 `lm_solve()`。],
)

返回值与 `LS_OneBreak` 一致：`coef`、`fitted`、`resi`、`SSE`。

== lm_solve

内部最小二乘求解器，对给定设计矩阵直接调用 `.lm.fit()`，避免构造并求逆
正规方程。`only.SSE = TRUE` 时仅返回残差平方和，用于断点搜索等只关心
拟合优度的场景。

#table(
  columns: (auto, 1fr),
  [参数], [说明],
  [`design`], [设计矩阵。],
  [`response`], [响应矩阵。],
  [`only.SSE`], [为 `TRUE` 时只返回 `SSE`。],
  [`...`], [保留参数，透传扩展。],
)

返回值：

+ `only.SSE = FALSE`：`coef`、`fitted`、`resi`、`SSE`；
+ `only.SSE = TRUE`：仅 `SSE`。
