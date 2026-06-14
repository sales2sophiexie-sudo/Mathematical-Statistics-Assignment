# Mathematical-Statistics-Assignment
========================================
项目名称：油气市场地缘政治事件冲击分析
========================================

1. 项目简介
-----------
本项目基于FRED原油价格、产量等日度数据（2025年6月–2026年6月），
采用三种互补方法评估地缘政治事件对WTI原油市场的冲击：

  方法A：Bai-Perron 多重结构断点检验（思路一）
  方法B：事件研究法 + 中断时间序列 + 冲击分离（思路二）

旨在识别市场结构性变化，量化事件影响，并对事件显著性进行排序。

2. 文件说明
-----------
(1) catch data.py
    - 功能：从FRED API获取原始数据
    - 输入：FRED API Key（已内置）
    - 输出：oil_analysis_data_2025_2026.csv
    - 使用方法：python catch data.py

(2) 1333.r（思路一：Bai-Perron断点检验）
    - 功能：对价格、收益率、产量序列进行真正的多重断点检验
    - 输出：
        * bai_perron_results.txt       （文本结果）
        * output_price_breaks.png      （价格断点图）
        * output_return_breaks.png     （收益率断点图）
        * output_production_breaks.png （产量断点图）
    - 使用方法：Rscript 1333.r

(3) 2.r（思路二：事件驱动分析）
    - 功能：
        * 事件研究法（CAR）
        * 中断时间序列（ITS，仅对E4/E8）
        * 分布滞后模型冲击分离
        * 综合显著性排序
    - 输出：
        * event_driven_results.txt           （详细结果文本）
        * output_car_barchart.png            （CAR条形图）
        * output_impact_decomposition.png    （冲击分解图）
        * output_composite_ranking.png       （综合得分排名图）
    - 使用方法：Rscript 2.r

3. 数据依赖
-----------
- FRED数据依赖如下series_id：
    POILWTIUSDM  : WTI原油价格（美元/桶）
    IPG211S      : 原油开采工业生产指数
    NASDAQQUSOI  : 原油进口价格指数（保留备用）

- R包依赖（思路一）：
    dplyr, tidyr, lubridate, readr, ggplot2, purrr

- R包依赖（思路二）：
    dplyr, tidyr, lubridate, readr, ggplot2, broom, sandwich, lmtest

- Python依赖：
    requests, pandas, datetime

4. 运行顺序（推荐）
-------------------
步骤1：python catch data.py        # 生成输入数据csv
步骤2：Rscript 1333.r              # 断点检验
步骤3：Rscript 2.r                 # 事件驱动分析

三个脚本可独立运行，但必须先执行步骤1生成数据文件。

5. 主要分析事件（思路二）
-------------------------
E1  2025-06-13  以空袭核设施        goldstein=-10
E2  2025-06-21  午夜之锤            goldstein=-10
E3  2025-06-23  停火协议            goldstein=+5
E4  2026-02-28  哈梅内伊身亡        goldstein=-10
E5  2026-03-01  伊朗报复            goldstein=-10
E6  2026-04-06  132次打击           goldstein=-10
E7  2026-06-07  导弹打击            goldstein=-10
E8  2026-06-10  封锁海峡            goldstein=-7

注：Goldstein分数基于CAMEO冲突编码（负=冲突升级，正=缓和）。

6. 关键参数说明（可自行调整）
------------------------------
- 事件研究法估计窗口：[-60, -11]
- 事件窗口：[-3, +7]
- ITS分析窗口：事件前后各60天
- 冲击分离半衰期：3天
- 分布滞后最大阶数：5
- Bai-Perron最小段长度：20–25天
- 断点-事件匹配容差：5天

7. 输出结果解读要点
--------------------
- CAR（累积异常收益）：>0为正冲击，<0为负冲击
- ITS：
    * Level Change：干预后是否立即跳变
    * Trend Change：长期趋势是否改变
- 冲击分离：每单位Goldstein对收益率的累计百分比影响
- 综合得分：结合CAR显著性与冲击贡献的排序指标

8. 注意事项
------------
- 数据为模拟演示数据（基于公开FRED代码结构），实际运行时需确保API有效。
- 部分ITS分析仅对E4/E8执行（因数据窗口限制）。
- 若时间序列自相关较强，建议增加Newey-West滞后阶数。
- 所有图片以PNG格式保存于当前工作目录。

9. 输出文件汇总
----------------
执行完整流程后，工作目录下将生成：
  oil_analysis_data_2025_2026.csv   （原始数据）
  bai_perron_results.txt
  output_price_breaks.png
  output_return_breaks.png
  output_production_breaks.png
  event_driven_results.txt
  output_car_barchart.png
  output_impact_decomposition.png
  output_composite_ranking.png

10. 联系与扩展建议
-------------------
- 可将ITS分析扩展至全部8个事件
- 可引入对照组（如其他商品指数）增强因果推断
- 可改用ARIMA模型处理残差自相关
- 可集成更多FRED变量（如库存、美元指数）

========================================
文档版本：1.0
适用脚本：catch data.py / 1333.r / 2.r
========================================
