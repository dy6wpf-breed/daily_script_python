import pandas as pd
import os
from datetime import datetime

# 定义文件路径
source_file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-凤仪场-只反馈DD和EE的填充基因型\【42】 KPSNY032024002C_KPS202502101-2_德康-凤仪种场_朱康平&罗毅&陈鹏宇&李肖晗_56个-20250514\KPSNY032024002C_朱康平、罗毅、陈鹏宇、李肖晗_德康-凤仪种场_56个_中芯一号PLUS_report_20250513\1.Sampleinfo\凤仪56个芯片位置对应表.xlsx"
# 读取Excel文件，跳过前面的行直到找到正确的表头
source_df = pd.read_excel(source_file_path, header=4)

# 定义需要的列顺序和名称
required_columns = [
    'sample_id', 'sample_name', 'sex', 'year', 'sire', 'dam', 'breed',
    '批次', '检测年份', '场', '芯片类型', '任务单号', '合同号', '实验室编号'
]

# 初始化一个新的DataFrame
extracted_data = pd.DataFrame()

# 添加缺失的列并用空字符串填充
for col in required_columns:
    extracted_data[col] = ''

# 提取存在的列
existing_columns = {
    'sample_id': 'Sample_ID',
    'sample_name': 'Sample_Name',
    '场': '客户单位',
    '任务单号': '任务单号',
    '合同号': '合同号',
    '实验室编号': '实验室编号'
}

for key, value in existing_columns.items():
    if value in source_df.columns:
        extracted_data[key] = source_df[value]

# 设置“芯片类型”为“Plus”
extracted_data['芯片类型'] = 'Plus'

# 设置“检测年份”为当前年份
current_year = datetime.now().year
extracted_data['检测年份'] = current_year

# 定义输出文件路径
output_folder = os.path.dirname(source_file_path)
output_file_name = '提取结果.xlsx'
output_file_path = os.path.join(output_folder, output_file_name)

# 将提取的数据写入新的Excel文件
extracted_data.to_excel(output_file_path, index=False)

print(f'提取完成，结果保存在：{output_file_path}')