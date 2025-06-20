import pandas as pd

# 读取td.csv文件
td_df = pd.read_csv('y:\\种质资源数据管理平台2024\\9-自贡德康\\td.csv')

# 读取Excel文件
excel_df = pd.read_excel('y:\\种质资源数据管理平台2024\\9-自贡德康\\1-样本基因型对应关系表模版-v1plus.xlsx')

# 匹配Sample_ID
matched_df = excel_df[excel_df['Sample_ID'].isin(td_df['Sample_ID'])]

# 将数据分成多份，每份最多9999条
chunk_size = 9999
chunks = [matched_df[i:i + chunk_size] for i in range(0, len(matched_df), chunk_size)]

# 保存分割后的数据到Excel，并提取Sample_ID到txt文件
for i, chunk in enumerate(chunks, 1):
    excel_output_path = f'y:\\种质资源数据管理平台2024\\9-自贡德康\\1-样本基因型对应关系表{i}.xlsx'
    txt_output_path = f'y:\\种质资源数据管理平台2024\\9-自贡德康\\1-样本基因型对应关系表{i}_Sample_IDs.txt'
    
    # 保存到Excel
    chunk.to_excel(excel_output_path, index=False)
    print(f"已保存Excel文件: {excel_output_path}")
    
    # 提取Sample_ID并保存到txt文件
    sample_ids = chunk['Sample_ID']
    sample_ids.to_csv(txt_output_path, index=False, header=False)
    print(f"已提取Sample_ID到文本文件: {txt_output_path}")

print("所有文件处理完成。")