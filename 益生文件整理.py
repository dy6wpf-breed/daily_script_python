import pandas as pd
import os
from datetime import datetime


# 定义文件路径，根据实际情况修改
data_path = r'Z:\猪中芯一号V1PLUS\【1】GS分析\【2】黑龙江益生-基因型数据TOP链要先转至GGP50K\【78】 KPSNY132023056C_KPS202501051_祝永华_黑龙江益生种猪繁育有限公司2506-2507_29个-20250116\KPSNY132023056C_祝永华_黑龙江益生种猪繁育有限公司2506-2507_29个_中芯一号PLUS_report_20250115\1.Sampleinfo'
result_path = r'Z:\猪中芯一号V1PLUS\【1】GS分析\【2】黑龙江益生-基因型数据TOP链要先转至GGP50K\【78】 KPSNY132023056C_KPS202501051_祝永华_黑龙江益生种猪繁育有限公司2506-2507_29个-20250116\KPSNY132023056C_祝永华_黑龙江益生种猪繁育有限公司2506-2507_29个_中芯一号PLUS_report_20250115\3.Analysis'


# 读取相关文件
excel_file_1 = pd.ExcelFile(f'{data_path}/2506-2507批次芯片位置对应表.xlsx')
df_1 = excel_file_1.parse('Sheet1', header=4)

excel_file_2 = pd.ExcelFile(f'{data_path}/系谱2506-07批次.XLSX')
df_2 = excel_file_2.parse('Sheet1')


# 获取 2506-2507批次2506-2507批次芯片位置对应表.xlsx 的创建时间
creation_time = os.path.getmtime(f'{data_path}/2506-2507批次芯片位置对应表.xlsx')
arrive_date = datetime.fromtimestamp(creation_time).strftime('%Y/%m/%d')


# 对于 1样本信息表-29_新.xlsx
df_sample = pd.DataFrame({
    'projectNo': df_1['合同号'],
    'taskNo': df_1['任务单号'],
    'enterpriseCode': 'SDYS4',
    'laboratoryNo': df_1['实验室编号'],
    'arriveDate': arrive_date,  # 使用文件创建时间
    'breed': 'YY',
    'sampleId': df_1['Sample_Name'],
    'populationInfo': '候选群'
})


# 对于 2 芯片位置对应表-29_新.xlsx
df_chip = pd.DataFrame({
    'sampleId': df_1['Sample_Name'],
    'laboratorySampleId': df_1['Sample_Name'],
    'chipPostionCode': df_1['Sample_ID'],
    'laboratoryNo': df_1['实验室编号'],
    'enterpriseCode': 'SDYS4',
    'populationInfo': '候选群',
    'breed': 'YY'
})


# 首先将df_2中的个体号与df_1中的Sample_Name进行匹配
df_2_matched = df_2[df_2['个体号'].isin(df_1['Sample_Name'])]

df_individual = pd.DataFrame({
    'sampleId': df_2_matched['个体号'],
    'breed': 'YY',
    'sex': df_2_matched['性别'].map({'公': 'M', '母': 'F'}),
    'strain': 'YY03',
    'birthDate': pd.to_datetime(df_2_matched['出生日期']).dt.strftime('%Y/%m/%d'),  # 修改日期格式
    'sire': df_2_matched['父亲'],
    'dam': df_2_matched['母亲']
})






# 计算样本数量，这里假设三个表的样本数量相同，取 df_1 的行数作为样本数量
sample_num = len(df_1)


# 保存为新文件，文件名加上样本数
df_sample.to_excel(f'{result_path}/1样本信息表-{sample_num}_新.xlsx', index=False)
df_chip.to_excel(f'{result_path}/2 芯片位置对应表-{sample_num}_新.xlsx', index=False)
df_individual.to_excel(f'{result_path}/3 个体信息表-{sample_num}_新.xlsx', index=False)