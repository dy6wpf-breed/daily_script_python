import pandas as pd

# 读取文本文件,只提取第一列
txt_file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-4个场数据\提取ID.txt"
txt_data = pd.read_csv(txt_file_path, sep="\t", usecols=[0], header=None, names=['first_column'])

# 读取Excel文件
excel_file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-4个场数据\德康四场ace文件全部-20250618.xlsx"
excel_data = pd.read_excel(excel_file_path)

# 进行匹配
merged_data = pd.merge(txt_data, excel_data, left_on='first_column', right_on='sample_name', how='left')

# 定义需要的表头顺序
new_column_order = ["first_column", "sample_id", "实验室编号", "sample_name", "任务单号", "sex", "year", "sire", "dam", "breed", "批次", "检测年份", "场", "芯片类型", "合同号", "系谱更正时间", "更改备注", "原始父", "原始母"]

# 调整列顺序
merged_data = merged_data[new_column_order]

# 输出结果
output_file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-4个场数据\德康四场错误ID结果.xlsx"
merged_data.to_excel(output_file_path, index=False)

print(f"匹配结果已保存到: {output_file_path}")
    