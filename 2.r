# ============================================================
# 思路二：事件驱动 – 外生事件显著性检验
# 方法：事件研究法(CAR) + 中断时间序列(ITS) + 冲击分离
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)
library(broom)
library(sandwich)
library(lmtest)

# ============================================================
# 1. 数据读取与预处理
# ============================================================

df_raw <- read_csv("oil_analysis_data_2025_2026.csv", show_col_types = FALSE)

df <- df_raw %>%
  pivot_wider(names_from = Variable, values_from = Value) %>%
  select(Date, 
         WTI_Price = `Global price of WTI Crude (美元/桶)`,
         Oil_Production = `Industrial Production: Mining: Crude Oil (指数)`) %>%
  mutate(Date = ymd(Date)) %>%
  arrange(Date)

# 填充完整日期
full_dates <- tibble(Date = seq(ymd("2025-06-01"), ymd("2026-06-11"), by = "day"))
df <- full_dates %>%
  left_join(df, by = "Date") %>%
  fill(WTI_Price, Oil_Production, .direction = "down") %>%
  drop_na() %>%
  mutate(WTI_Return = (WTI_Price - lag(WTI_Price)) / lag(WTI_Price) * 100) %>%
  drop_na()

cat("数据范围:", as.character(min(df$Date)), "至", as.character(max(df$Date)))
cat("\n总天数:", nrow(df), "\n")

# ============================================================
# 2. 定义事件（Goldstein分数基于CAMEO表）
# ============================================================

events <- data.frame(
  label = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8"),
  date = as.Date(c("2025-06-13", "2025-06-21", "2025-06-23",
                   "2026-02-28", "2026-03-01", "2026-04-06",
                   "2026-06-07", "2026-06-10")),
  desc = c("以空袭核设施", "午夜之锤", "停火协议",
           "哈梅内伊身亡", "伊朗报复", "132次打击",
           "导弹打击", "封锁海峡"),
  goldstein = c(-10.0, -10.0, 5.0, -10.0, -10.0, -10.0, -10.0, -7.0)
)

cat("\n已定义事件:\n")
print(events[, c("label", "date", "goldstein")])

# ============================================================
# 3. 事件研究法 - 计算累积异常收益 (CAR)（修正版）
# ============================================================

calculate_car <- function(returns, event_idx, est_window = c(-60, -11), 
                          event_window = c(-3, 7)) {
  
  # 初始化返回值
  result <- list(
    car = NA_real_,
    t_stat = NA_real_,
    p_value = NA_real_,
    is_significant = FALSE,
    sig_level = "n.s."
  )
  
  # 检查输入有效性
  if (is.na(event_idx) || event_idx < 1 || event_idx > length(returns)) {
    return(result)
  }
  
  est_start <- event_idx + est_window[1]
  est_end <- event_idx + est_window[2]
  
  # 检查估计窗口边界
  if (est_start < 1 || est_end > length(returns) || est_start >= est_end) {
    return(result)
  }
  
  # 正常收益（估计窗口均值）
  mu <- mean(returns[est_start:est_end], na.rm = TRUE)
  sigma <- sd(returns[est_start:est_end], na.rm = TRUE)
  
  if (is.na(sigma) || sigma == 0) {
    return(result)
  }
  
  # 异常收益
  ar <- c()
  for (t in event_window[1]:event_window[2]) {
    idx <- event_idx + t
    if (idx >= 1 && idx <= length(returns)) {
      ar <- c(ar, returns[idx] - mu)
    }
  }
  
  if (length(ar) == 0) {
    return(result)
  }
  
  # 累积异常收益
  car <- sum(ar)
  std_car <- sqrt(length(ar)) * sigma
  
  if (is.na(std_car) || std_car == 0) {
    return(result)
  }
  
  # t检验
  t_stat <- car / std_car
  p_value <- 2 * (1 - pt(abs(t_stat), df = length(ar) - 1))
  
  result$car <- car
  result$t_stat <- t_stat
  result$p_value <- p_value
  result$is_significant <- p_value < 0.05
  result$sig_level <- ifelse(p_value < 0.01, "***", 
                             ifelse(p_value < 0.05, "**", 
                                    ifelse(p_value < 0.1, "*", "n.s.")))
  
  return(result)
}

# 对每个事件计算CAR（修正版：逐行添加，避免data.frame长度问题）
car_results <- data.frame()

for (i in 1:nrow(events)) {
  event_date <- events$date[i]
  event_idx <- which(df$Date == event_date)
  
  cat("处理事件:", events$label[i], "日期:", as.character(event_date), "\n")
  
  if (length(event_idx) == 0) {
    cat("  警告: 事件日期不在数据范围内，跳过\n")
    next
  }
  
  # 油价CAR
  car_price <- calculate_car(df$WTI_Return, event_idx)
  # 产量CAR
  car_prod <- calculate_car(df$Oil_Production, event_idx)
  
  # 使用rbind逐行添加
  new_row <- data.frame(
    Event = events$label[i],
    Date = as.character(event_date),
    Desc = events$desc[i],
    Goldstein = events$goldstein[i],
    CAR_Price = car_price$car,
    Sig_Price = car_price$sig_level,
    P_Price = car_price$p_value,
    CAR_Prod = car_prod$car,
    Sig_Prod = car_prod$sig_level,
    P_Prod = car_prod$p_value,
    stringsAsFactors = FALSE
  )
  
  car_results <- bind_rows(car_results, new_row)
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("事件研究法 (CAR) 结果\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
print(car_results[, c("Event", "Desc", "CAR_Price", "Sig_Price", "CAR_Prod", "Sig_Prod")])

# ============================================================
# 4. 中断时间序列分析 (ITS) - 检验水平变化和趋势变化
# ============================================================

its_analysis <- function(df_sub, event_date, variable) {
  # 初始化返回值
  result <- list(
    level_change = NA_real_,
    level_pval = NA_real_,
    level_sig = "n.s.",
    trend_change = NA_real_,
    trend_pval = NA_real_,
    trend_sig = "n.s."
  )
  
  if (nrow(df_sub) < 20) {
    return(result)
  }
  
  # 创建时间变量
  df_sub <- df_sub %>%
    mutate(
      Time = row_number(),
      Post = ifelse(Date >= event_date, 1, 0),
      Time_Post = Time * Post
    )
  
  # 回归模型
  formula <- as.formula(paste(variable, "~ Time + Post + Time_Post"))
  model <- tryCatch({
    lm(formula, data = df_sub)
  }, error = function(e) NULL)
  
  if (is.null(model)) {
    return(result)
  }
  
  # Newey-West HAC标准误
  vcov_hac <- tryCatch({
    NeweyWest(model, lag = 5, prewhite = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(vcov_hac)) {
    return(result)
  }
  
  coeftest_result <- coeftest(model, vcov = vcov_hac)
  
  # 提取结果
  if ("Post" %in% rownames(coeftest_result)) {
    result$level_change <- coeftest_result["Post", "Estimate"]
    result$level_pval <- coeftest_result["Post", "Pr(>|t|)"]
    result$level_sig <- ifelse(result$level_pval < 0.01, "***",
                               ifelse(result$level_pval < 0.05, "**",
                                      ifelse(result$level_pval < 0.1, "*", "n.s.")))
  }
  
  if ("Time_Post" %in% rownames(coeftest_result)) {
    result$trend_change <- coeftest_result["Time_Post", "Estimate"]
    result$trend_pval <- coeftest_result["Time_Post", "Pr(>|t|)"]
    result$trend_sig <- ifelse(result$trend_pval < 0.01, "***",
                               ifelse(result$trend_pval < 0.05, "**",
                                      ifelse(result$trend_pval < 0.1, "*", "n.s.")))
  }
  
  return(result)
}

# 只对可能产生结构变化的事件进行ITS（E4和E8）
its_results <- data.frame()

for (event_label in c("E4", "E8")) {
  event_info <- events[events$label == event_label, ]
  event_date <- event_info$date
  window_days <- 60
  
  # 截取事件前后60天的窗口
  df_sub <- df %>%
    filter(Date >= event_date - days(window_days),
           Date <= event_date + days(window_days))
  
  cat("\n处理ITS事件:", event_label, "窗口大小:", nrow(df_sub), "天\n")
  
  if (nrow(df_sub) < 20) next
  
  # 油价ITS
  its_price <- its_analysis(df_sub, event_date, "WTI_Price")
  # 产量ITS
  its_prod <- its_analysis(df_sub, event_date, "Oil_Production")
  
  new_row <- data.frame(
    Event = event_label,
    Desc = event_info$desc,
    Price_Level_Change = its_price$level_change,
    Price_Level_Sig = its_price$level_sig,
    Price_Trend_Change = its_price$trend_change,
    Price_Trend_Sig = its_price$trend_sig,
    Prod_Level_Change = its_prod$level_change,
    Prod_Level_Sig = its_prod$level_sig,
    stringsAsFactors = FALSE
  )
  
  its_results <- bind_rows(its_results, new_row)
}

cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
cat("中断时间序列分析 (ITS) 结果 - E4和E8\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
if (nrow(its_results) > 0) {
  print(its_results)
} else {
  cat("未生成ITS结果\n")
}

# ============================================================
# 5. 冲击分离（分布滞后模型 + Goldstein）
# ============================================================

cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
cat("冲击分离（分布滞后模型）\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 构建每日Goldstein强度序列（指数衰减叠加）
half_life <- 3
lambda_decay <- log(2) / half_life

df$Goldstein_Intensity <- 0

for (i in 1:nrow(events)) {
  event_date <- events$date[i]
  goldstein <- events$goldstein[i]
  
  # 计算衰减权重
  days_after <- as.numeric(difftime(df$Date, event_date, units = "days"))
  decay <- ifelse(days_after >= 0, exp(-lambda_decay * days_after), 0)
  df$Goldstein_Intensity <- df$Goldstein_Intensity + goldstein * decay
}

# 分布滞后模型
max_lag <- 5

# 创建滞后变量
for (lag in 0:max_lag) {
  df[[paste0("GS_lag", lag)]] <- lag(df$Goldstein_Intensity, lag)
}

# 删除缺失值
df_lag <- df[complete.cases(df[, paste0("GS_lag", 0:max_lag)]), ]

if (nrow(df_lag) > max_lag * 2) {
  # 油价模型
  formula_price <- as.formula(paste("WTI_Return ~", paste(paste0("GS_lag", 0:max_lag), collapse = " + ")))
  model_price <- lm(formula_price, data = df_lag)
  
  vcov_hac_price <- tryCatch({
    NeweyWest(model_price, lag = max_lag, prewhite = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(vcov_hac_price)) {
    coeftest_price <- coeftest(model_price, vcov = vcov_hac_price)
    cum_effect_price <- sum(coeftest_price[2:(max_lag+2), "Estimate"])
  } else {
    cum_effect_price <- sum(coef(model_price)[2:(max_lag+2)])
  }
  
  # 产量模型
  formula_prod <- as.formula(paste("Oil_Production ~", paste(paste0("GS_lag", 0:max_lag), collapse = " + ")))
  model_prod <- lm(formula_prod, data = df_lag)
  
  vcov_hac_prod <- tryCatch({
    NeweyWest(model_prod, lag = max_lag, prewhite = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(vcov_hac_prod)) {
    coeftest_prod <- coeftest(model_prod, vcov = vcov_hac_prod)
    cum_effect_prod <- sum(coeftest_prod[2:(max_lag+2), "Estimate"])
  } else {
    cum_effect_prod <- sum(coef(model_prod)[2:(max_lag+2)])
  }
  
  cat("\n累计冲击效应（每单位Goldstein）:")
  cat(paste("\n  油价:", round(cum_effect_price, 4), "百分点"))
  cat(paste("\n  产量:", round(cum_effect_prod, 4), "百分点"))
  
  # 计算每个事件的边际贡献
  contrib_results <- events %>%
    mutate(
      Price_Contribution = goldstein * cum_effect_price,
      Prod_Contribution = goldstein * cum_effect_prod
    ) %>%
    arrange(desc(Price_Contribution))
  
  cat("\n\n每个事件的边际贡献:\n")
  print(contrib_results[, c("label", "desc", "goldstein", "Price_Contribution", "Prod_Contribution")])
  
} else {
  cat("\n样本量不足，无法进行分布滞后模型分析\n")
  cum_effect_price <- NA
  cum_effect_prod <- NA
  contrib_results <- events %>%
    mutate(
      Price_Contribution = NA,
      Prod_Contribution = NA
    )
}

# ============================================================
# 6. 综合显著性排序
# ============================================================

if (exists("cum_effect_price") && !is.na(cum_effect_price) && nrow(car_results) > 0) {
  
  final_results <- car_results %>%
    left_join(contrib_results[, c("label", "Price_Contribution")], by = c("Event" = "label")) %>%
    mutate(
      # 显著性得分：***=3, **=2, *=1, n.s.=0
      Sig_Score = case_when(
        Sig_Price == "***" ~ 3,
        Sig_Price == "**" ~ 2,
        Sig_Price == "*" ~ 1,
        TRUE ~ 0
      ),
      # 综合得分
      Composite_Score = Sig_Score + (Price_Contribution - min(Price_Contribution, na.rm = TRUE)) / 
                        (max(Price_Contribution, na.rm = TRUE) - min(Price_Contribution, na.rm = TRUE)) * 2
    ) %>%
    arrange(desc(Composite_Score))
  
  cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("综合显著性排序\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  print(final_results[, c("Event", "Desc", "CAR_Price", "Sig_Price", "Price_Contribution", "Composite_Score")])
  
} else {
  cat("\n无法计算综合排序（数据不足）\n")
}

# ============================================================
# 7. 可视化（三张独立图片）
# ============================================================

# 图1: CAR条形图
if (nrow(car_results) > 0) {
  p1 <- ggplot(car_results, aes(x = reorder(Event, -CAR_Price), y = CAR_Price, fill = Sig_Price)) +
    geom_col() +
    geom_hline(yintercept = 0, color = "black") +
    scale_fill_manual(values = c("***" = "darkred", "**" = "red", "*" = "salmon", "n.s." = "gray")) +
    labs(title = "Cumulative Abnormal Return (CAR) by Event",
         x = "Event", y = "CAR (%)", fill = "Significance") +
    theme_minimal()
  ggsave("output_car_barchart.png", p1, width = 10, height = 6, dpi = 150)
  cat("\n✓ 图1已保存: output_car_barchart.png\n")
}

# 图2: 冲击分离贡献
if (exists("contrib_results") && nrow(contrib_results) > 0 && !all(is.na(contrib_results$Price_Contribution))) {
  p2 <- ggplot(contrib_results, aes(x = reorder(label, -Price_Contribution), 
                                     y = Price_Contribution, fill = label)) +
    geom_col(show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "black") +
    labs(title = "Decomposed Impact per Event (Goldstein × Cumulative Effect)",
         x = "Event", y = "Contribution (%)") +
    theme_minimal()
  ggsave("output_impact_decomposition.png", p2, width = 10, height = 6, dpi = 150)
  cat("✓ 图2已保存: output_impact_decomposition.png\n")
}

# 图3: 综合得分排名
if (exists("final_results") && nrow(final_results) > 0) {
  p3 <- ggplot(final_results, aes(x = reorder(Event, Composite_Score), 
                                   y = Composite_Score, fill = Sig_Price)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("***" = "darkred", "**" = "red", "*" = "salmon", "n.s." = "gray")) +
    labs(title = "Composite Score Ranking (Higher = More Significant)",
         x = "Event", y = "Composite Score", fill = "Significance") +
    theme_minimal()
  ggsave("output_composite_ranking.png", p3, width = 10, height = 6, dpi = 150)
  cat("✓ 图3已保存: output_composite_ranking.png\n")
}

# ============================================================
# 8. 保存结果到文本文件
# ============================================================

sink("event_driven_results.txt")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("思路二：事件驱动分析结果\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("【1. 事件研究法 (CAR) 结果】\n")
if (nrow(car_results) > 0) {
  print(car_results[, c("Event", "Desc", "CAR_Price", "Sig_Price", "CAR_Prod", "Sig_Prod")])
} else {
  cat("无CAR结果\n")
}

cat("\n\n【2. 中断时间序列分析 (ITS) 结果】\n")
if (nrow(its_results) > 0) {
  print(its_results)
} else {
  cat("无ITS结果\n")
}

cat("\n\n【3. 冲击分离结果】\n")
if (exists("cum_effect_price") && !is.na(cum_effect_price)) {
  cat(paste("累计冲击效应（油价）:", round(cum_effect_price, 4), "%/每单位Goldstein\n"))
  cat(paste("累计冲击效应（产量）:", round(cum_effect_prod, 4), "%/每单位Goldstein\n\n"))
  print(contrib_results[, c("label", "desc", "goldstein", "Price_Contribution", "Prod_Contribution")])
} else {
  cat("无冲击分离结果\n")
}

cat("\n\n【4. 综合显著性排序】\n")
if (exists("final_results") && nrow(final_results) > 0) {
  print(final_results[, c("Event", "Desc", "CAR_Price", "Sig_Price", "Price_Contribution", "Composite_Score")])
} else {
  cat("无综合排序结果\n")
}

sink()

cat("\n\n✅ 思路二分析完成！")
cat("\n生成的文件:")
cat("\n  结果文本: event_driven_results.txt")
cat("\n  图1(CAR): output_car_barchart.png")
cat("\n  图2(冲击分离): output_impact_decomposition.png")
cat("\n  图3(综合得分): output_composite_ranking.png\n")