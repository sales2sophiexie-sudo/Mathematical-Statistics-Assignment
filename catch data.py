import requests
import pandas as pd
from datetime import datetime

# 你的 FRED API 密钥
api_key = "528d1e6abeb0535a6cb42dbc1d278aa7"

# 定义指标代码 (Series ID) 和中文名称映射
# 经过核实，以下是符合你需求的标准FRED代码：
series_map = {
    "POILWTIUSDM": "Global price of WTI Crude (美元/桶)",
    "IPG211S": "Industrial Production: Mining: Crude Oil (指数)",
    "NASDAQQUSOI": "Import Price Index (End Use): Crude Oil (指数)"
}

# 时间范围：2025.6 - 2026.6
start_date = "2025-06-01"
end_date = "2026-06-30"

# 存储结果的字典
all_data = {}

for series_id, name in series_map.items():
    url = "https://api.stlouisfed.org/fred/series/observations"
    params = {
        "series_id": series_id,
        "api_key": api_key,
        "file_type": "json",
        "observation_start": start_date,
        "observation_end": end_date,
        "units": "lin"  # 线性刻度，获取原始数据
    }

    response = requests.get(url, params=params)

    if response.status_code == 200:
        data = response.json()
        # 提取日期和数值
        observations = []
        for obs in data['observations']:
            date_str = obs['date']
            value = obs['value']
            # 过滤掉缺失值 "."
            if value != ".":
                observations.append({"Date": date_str, "Variable": name, "Value": float(value)})

        all_data[series_id] = observations
        print(f"成功获取数据: {name}")
    else:
        print(f"获取数据失败: {name}, 状态码: {response.status_code}")

# 将数据合并并保存为 CSV (方便导入统计软件)
combined_list = []
for series_id, obs_list in all_data.items():
    combined_list.extend(obs_list)

df = pd.DataFrame(combined_list)
df.sort_values(['Date', 'Variable'], inplace=True)
print("\n数据预览:")
print(df.head(10))

# 保存到文件
output_file = "oil_analysis_data_2025_2026.csv"
import os
print("当前文件将保存在:", os.getcwd())
df.to_csv(output_file, index=False)
print(f"\n数据已保存至: {output_file}")