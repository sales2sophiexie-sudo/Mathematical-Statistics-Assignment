# ============================================================
# 思路一：真正的 Bai-Perron 多重结构断点检验
# 实现：多断点搜索 + UDmax检验 + 序贯SupF检验
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)
library(purrr)

# 1. 数据读取与预处理
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

n <- nrow(df)
print(paste("样本长度:", n))

# ============================================================
# 2. 真正的 Bai-Perron 多重断点检验函数
# ============================================================

# 2.1 动态规划搜索最优多断点位置
dynamic_programming_breaks <- function(y, min_seg, max_breaks) {
  N <- length(y)
  
  # 初始化DP表: dp[i][k] = 前i个点、k个断点的最小RSS
  dp <- matrix(Inf, nrow = N + 1, ncol = max_breaks + 1)
  split_point <- matrix(0, nrow = N + 1, ncol = max_breaks + 1)
  
  # 0个断点的情况
  for (i in min_seg:N) {
    segment <- y[1:i]
    mean_seg <- mean(segment)
    dp[i, 1] <- sum((segment - mean_seg)^2)
  }
  
  # 动态规划递推
  for (k in 2:(max_breaks + 1)) {
    for (i in (k * min_seg):N) {
      for (j in ((k-1) * min_seg):(i - min_seg)) {
        segment <- y[(j+1):i]
        mean_seg <- mean(segment)
        rss_last <- sum((segment - mean_seg)^2)
        total_rss <- dp[j, k-1] + rss_last
        if (total_rss < dp[i, k]) {
          dp[i, k] <- total_rss
          split_point[i, k] <- j
        }
      }
    }
  }
  
  # 回溯获取最优断点位置
  get_breaks <- function(k) {
    if (k == 0) return(integer(0))
    breaks <- integer(k)
    current <- N
    for (i in k:1) {
      breaks[i] <- split_point[current, i + 1]
      current <- breaks[i]
    }
    return(breaks)
  }
  
  return(list(dp = dp, split_point = split_point, get_breaks = get_breaks))
}

# 2.2 计算给定断点位置的RSS
compute_rss <- function(y, break_positions) {
  N <- length(y)
  if (length(break_positions) == 0) {
    mean_all <- mean(y)
    return(sum((y - mean_all)^2))
  }
  
  break_positions <- sort(break_positions)
  segments <- c(0, break_positions, N)
  rss_total <- 0
  
  for (i in 1:(length(segments) - 1)) {
    start <- segments[i] + 1
    end <- segments[i + 1]
    segment <- y[start:end]
    mean_seg <- mean(segment)
    rss_total <- rss_total + sum((segment - mean_seg)^2)
  }
  return(rss_total)
}

# 2.3 UDmax检验（存在性检验）
udmax_test <- function(y, min_seg) {
  N <- length(y)
  overall_mean <- mean(y)
  rss_null <- sum((y - overall_mean)^2)
  
  f_stats <- numeric()
  candidate_pos <- integer()
  
  for (t in min_seg:(N - min_seg)) {
    seg1 <- y[1:t]
    seg2 <- y[(t+1):N]
    mean1 <- mean(seg1)
    mean2 <- mean(seg2)
    rss_alt <- sum((seg1 - mean1)^2) + sum((seg2 - mean2)^2)
    f_stat <- ((rss_null - rss_alt) / 1) / (rss_alt / (N - 2))
    f_stats <- c(f_stats, f_stat)
    candidate_pos <- c(candidate_pos, t)
  }
  
  udmax <- max(f_stats)
  best_pos <- candidate_pos[which.max(f_stats)]
  
  # Bai-Perron 5%临界值（近似）
  critical_5pct <- 8.85
  
  return(list(
    udmax = udmax,
    best_pos = best_pos,
    is_significant = udmax > critical_5pct,
    f_stats = f_stats,
    candidate_pos = candidate_pos
  ))
}

# 2.4 序贯检验（确定最优断点数量）
sequential_test <- function(y, min_seg, max_breaks = 5) {
  N <- length(y)
  results <- list()
  
  # 首先运行UDmax检验
  udmax_res <- udmax_test(y, min_seg)
  if (!udmax_res$is_significant) {
    return(list(optimal_breaks = 0, break_positions = integer(0), 
                is_sig = FALSE, f_stats = list()))
  }
  
  current_breaks <- integer(0)
  current_rss <- compute_rss(y, current_breaks)
  
  for (m in 0:(max_breaks - 1)) {
    # 尝试增加一个断点
    candidate_rss <- numeric()
    candidate_positions <- list()
    
    # 扫描每个可能的新断点位置
    for (t in (min_seg):(N - min_seg)) {
      # 确保新断点不与现有断点冲突
      if (any(abs(t - current_breaks) < min_seg)) next
      
      # 构造新的断点集合
      new_breaks <- sort(c(current_breaks, t))
      rss_alt <- compute_rss(y, new_breaks)
      candidate_rss <- c(candidate_rss, rss_alt)
      candidate_positions <- c(candidate_positions, list(new_breaks))
    }
    
    if (length(candidate_rss) == 0) break
    
    # 找到最佳新断点
    best_idx <- which.min(candidate_rss)
    best_new_breaks <- candidate_positions[[best_idx]]
    rss_alt <- candidate_rss[best_idx]
    
    # 计算F统计量（SupF(m+1|m)）
    df1 <- 2
    df2 <- N - 2 * (m + 1) - 1
    f_stat <- ((current_rss - rss_alt) / df1) / (rss_alt / df2)
    
    results[[paste0("m=", m)]] <- list(
      current_breaks = current_breaks,
      proposed_breaks = best_new_breaks,
      f_stat = f_stat,
      rss_null = current_rss,
      rss_alt = rss_alt
    )
    
    critical_5pct <- 8.85
    
    if (f_stat > critical_5pct) {
      # 显著，接受新断点
      current_breaks <- best_new_breaks
      current_rss <- rss_alt
    } else {
      # 不显著，停止
      return(list(optimal_breaks = m, break_positions = current_breaks,
                  is_sig = TRUE, f_stats = results))
    }
  }
  
  return(list(optimal_breaks = max_breaks, break_positions = current_breaks,
              is_sig = TRUE, f_stats = results))
}

# ============================================================
# 3. 执行真正的 Bai-Perron 检验
# ============================================================

min_seg_price <- 25
min_seg_return <- 20
min_seg_prod <- 25
max_breaks <- 3

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Bai-Perron 多重结构断点检验结果\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 3.1 价格水平
cat("\n【1. WTI 价格水平】\n")
udmax_price <- udmax_test(df$WTI_Price, min_seg_price)
cat(paste("  UDmax统计量:", round(udmax_price$udmax, 2)))
cat(paste("\n  是否有断点:", ifelse(udmax_price$is_significant, "是", "否")))

if (udmax_price$is_significant) {
  seq_price <- sequential_test(df$WTI_Price, min_seg_price, max_breaks)
  cat(paste("\n  最优断点数量:", seq_price$optimal_breaks))
  if (length(seq_price$break_positions) > 0) {
    break_dates_price <- df$Date[seq_price$break_positions]
    cat(paste("\n  断点日期:", paste(as.character(break_dates_price), collapse = ", ")))
  } else {
    break_dates_price <- as.Date(character(0))
  }
} else {
  seq_price <- list(optimal_breaks = 0, break_positions = integer(0))
  break_dates_price <- as.Date(character(0))
}

# 3.2 收益率
cat("\n\n【2. WTI 收益率】\n")
udmax_return <- udmax_test(df$WTI_Return, min_seg_return)
cat(paste("  UDmax统计量:", round(udmax_return$udmax, 2)))
cat(paste("\n  是否有断点:", ifelse(udmax_return$is_significant, "是", "否")))

if (udmax_return$is_significant) {
  seq_return <- sequential_test(df$WTI_Return, min_seg_return, max_breaks)
  cat(paste("\n  最优断点数量:", seq_return$optimal_breaks))
  if (length(seq_return$break_positions) > 0) {
    break_dates_return <- df$Date[seq_return$break_positions]
    cat(paste("\n  断点日期:", paste(as.character(break_dates_return), collapse = ", ")))
  } else {
    break_dates_return <- as.Date(character(0))
  }
} else {
  seq_return <- list(optimal_breaks = 0, break_positions = integer(0))
  break_dates_return <- as.Date(character(0))
}

# 3.3 工业产量
cat("\n\n【3. 工业产量指数】\n")
udmax_prod <- udmax_test(df$Oil_Production, min_seg_prod)
cat(paste("  UDmax统计量:", round(udmax_prod$udmax, 2)))
cat(paste("\n  是否有断点:", ifelse(udmax_prod$is_significant, "是", "否")))

if (udmax_prod$is_significant) {
  seq_prod <- sequential_test(df$Oil_Production, min_seg_prod, max_breaks)
  cat(paste("\n  最优断点数量:", seq_prod$optimal_breaks))
  if (length(seq_prod$break_positions) > 0) {
    break_dates_prod <- df$Date[seq_prod$break_positions]
    cat(paste("\n  断点日期:", paste(as.character(break_dates_prod), collapse = ", ")))
  } else {
    break_dates_prod <- as.Date(character(0))
  }
} else {
  seq_prod <- list(optimal_breaks = 0, break_positions = integer(0))
  break_dates_prod <- as.Date(character(0))
}

# ============================================================
# 4. 已知事件匹配分析
# ============================================================

events <- data.frame(
  label = c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8"),
  date = as.Date(c("2025-06-13", "2025-06-21", "2025-06-23",
                   "2026-02-28", "2026-03-01", "2026-04-06",
                   "2026-06-07", "2026-06-10")),
  desc = c("以空袭核设施", "午夜之锤", "停火协议",
           "哈梅内伊身亡", "伊朗报复", "132次打击",
           "导弹打击", "封锁海峡")
)

TOLERANCE <- 5
cat("\n\n", paste(rep("=", 70), collapse = ""), "\n")
cat(paste("断点与已知事件匹配分析（容差", TOLERANCE, "天）\n"))
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\n| 事件 | 日期 | 描述 | 价格断点 | 收益率断点 | 产量断点 |\n")
cat("|------|------|------|----------|------------|----------|\n")

for (i in 1:nrow(events)) {
  event_date <- events$date[i]
  
  price_match <- any(abs(as.numeric(event_date - break_dates_price)) <= TOLERANCE)
  return_match <- any(abs(as.numeric(event_date - break_dates_return)) <= TOLERANCE)
  prod_match <- any(abs(as.numeric(event_date - break_dates_prod)) <= TOLERANCE)
  
  price_str <- ifelse(price_match, "✓", "✗")
  return_str <- ifelse(return_match, "✓", "✗")
  prod_str <- ifelse(prod_match, "✓", "✗")
  
  cat(paste("|", events$label[i], "|", as.character(event_date), 
            "|", events$desc[i], "|", price_str, "|", return_str, "|", prod_str, "|\n"))
}

# ============================================================
# 5. 可视化（三张独立图片）
# ============================================================

# 图1: 价格水平
p1 <- ggplot(df, aes(x = Date, y = WTI_Price)) +
  geom_line(color = "black", linewidth = 1) +
  { if(length(break_dates_price) > 0) 
      geom_vline(xintercept = as.numeric(break_dates_price), 
                 color = "red", linetype = "dashed", linewidth = 1) } +
  geom_vline(xintercept = as.numeric(events$date), 
             color = "purple", alpha = 0.6, linewidth = 0.8) +
  labs(title = paste0("WTI Price - Bai-Perron Breaks (Red) vs Events (Purple) | UDmax=", round(udmax_price$udmax, 2)),
       x = "Date", y = "Price (USD/barrel)") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("output_price_breaks.png", p1, width = 14, height = 6, dpi = 150)
cat("\n✓ 图1已保存: output_price_breaks.png")

# 图2: 收益率
p2 <- ggplot(df, aes(x = Date, y = WTI_Return)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  { if(length(break_dates_return) > 0) 
      geom_vline(xintercept = as.numeric(break_dates_return), 
                 color = "red", linetype = "dashed", linewidth = 1) } +
  geom_vline(xintercept = as.numeric(events$date), 
             color = "purple", alpha = 0.6, linewidth = 0.8) +
  labs(title = paste0("WTI Returns - Bai-Perron Breaks (Red) vs Events (Purple) | UDmax=", round(udmax_return$udmax, 2)),
       x = "Date", y = "Return (%)") +
  theme_minimal()
ggsave("output_return_breaks.png", p2, width = 14, height = 6, dpi = 150)
cat("\n✓ 图2已保存: output_return_breaks.png")

# 图3: 产量
p3 <- ggplot(df, aes(x = Date, y = Oil_Production)) +
  geom_line(color = "green", linewidth = 1) +
  { if(length(break_dates_prod) > 0) 
      geom_vline(xintercept = as.numeric(break_dates_prod), 
                 color = "red", linetype = "dashed", linewidth = 1) } +
  geom_vline(xintercept = as.numeric(events$date), 
             color = "purple", alpha = 0.6, linewidth = 0.8) +
  labs(title = paste0("Oil Production - Bai-Perron Breaks (Red) vs Events (Purple) | UDmax=", round(udmax_prod$udmax, 2)),
       x = "Date", y = "Production Index") +
  theme_minimal()
ggsave("output_production_breaks.png", p3, width = 14, height = 6, dpi = 150)
cat("\n✓ 图3已保存: output_production_breaks.png")

# ============================================================
# 6. 保存结果到文本文件
# ============================================================

sink("bai_perron_results.txt")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("真正的 Bai-Perron 多重结构断点检验结果\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("【1. WTI 价格水平】\n")
cat(paste("  UDmax统计量:", round(udmax_price$udmax, 2), "\n"))
cat(paste("  是否存在断点:", ifelse(udmax_price$is_significant, "是", "否"), "\n"))
cat(paste("  最优断点数量:", seq_price$optimal_breaks, "\n"))
if (length(break_dates_price) > 0) {
  cat(paste("  断点日期:", paste(as.character(break_dates_price), collapse = ", "), "\n"))
}

cat("\n【2. WTI 收益率】\n")
cat(paste("  UDmax统计量:", round(udmax_return$udmax, 2), "\n"))
cat(paste("  是否存在断点:", ifelse(udmax_return$is_significant, "是", "否"), "\n"))
cat(paste("  最优断点数量:", seq_return$optimal_breaks, "\n"))
if (length(break_dates_return) > 0) {
  cat(paste("  断点日期:", paste(as.character(break_dates_return), collapse = ", "), "\n"))
}

cat("\n【3. 工业产量指数】\n")
cat(paste("  UDmax统计量:", round(udmax_prod$udmax, 2), "\n"))
cat(paste("  是否存在断点:", ifelse(udmax_prod$is_significant, "是", "否"), "\n"))
cat(paste("  最优断点数量:", seq_prod$optimal_breaks, "\n"))
if (length(break_dates_prod) > 0) {
  cat(paste("  断点日期:", paste(as.character(break_dates_prod), collapse = ", "), "\n"))
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("断点与已知事件匹配分析（容差5天）\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

for (i in 1:nrow(events)) {
  event_date <- events$date[i]
  price_match <- any(abs(as.numeric(event_date - break_dates_price)) <= TOLERANCE)
  return_match <- any(abs(as.numeric(event_date - break_dates_return)) <= TOLERANCE)
  prod_match <- any(abs(as.numeric(event_date - break_dates_prod)) <= TOLERANCE)
  
  cat(paste(events$label[i], events$desc[i], "(", as.character(event_date), "): "))
  cat(paste("价格:", ifelse(price_match, "匹配", "不匹配")))
  cat(paste("; 收益率:", ifelse(return_match, "匹配", "不匹配")))
  cat(paste("; 产量:", ifelse(prod_match, "匹配", "不匹配"), "\n"))
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("分析完成！\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

sink()

cat("\n\n✅ 分析完成！")
cat("\n生成的文件:")
cat("\n  结果文本: bai_perron_results.txt")
cat("\n  图1(价格): output_price_breaks.png")
cat("\n  图2(收益率): output_return_breaks.png")
cat("\n  图3(产量): output_production_breaks.png")