import pandas as pd
import os

# 定义文件路径
excel_file = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-葛溪种场\【9】KPSNY032025019C_KPS202506054_德康-葛溪公猪站_朱康平&罗毅&陈鹏宇&李肖晗_29个-20250616\KPSNY032025019C_朱康平、罗毅、陈鹏宇、李肖晗_德康-葛溪公猪站_29个_中芯一号PLUS_report_20250616\1.Sampleinfo\葛溪29个芯片位置对应表.xlsx"
csv_file = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-葛溪种场\【9】KPSNY032025019C_KPS202506054_德康-葛溪公猪站_朱康平&罗毅&陈鹏宇&李肖晗_29个-20250616\download\2系谱纠偏结果\ped_correction.csv"
output_folder = r'Z:\猪中芯一号V1PLUS\【1】GS分析\【1】德康\德康-葛溪种场\【9】KPSNY032025019C_KPS202506054_德康-葛溪公猪站_朱康平&罗毅&陈鹏宇&李肖晗_29个-20250616\download\2系谱纠偏结果'
ped_error_file = os.path.join(output_folder, 'ped_error.csv')

try:
    # 从 Excel 文件的第 5 行开始读取数据
    df_excel = pd.read_excel(excel_file, header=4)

    # 读取 CSV 文件
    df_csv = pd.read_csv(csv_file, encoding='gbk')

    # 合并数据
    merged_df = pd.merge(df_excel[['Sample_ID', 'Sample_Name']], df_csv[['ID', 'Fa_ID', 'Mo_ID', 'Fa_infor', 'Mo_infor', 'ori_Fa_ID', 'ori_Mo_ID']], left_on='Sample_Name', right_on='ID', how='inner')

    # 删除 ID 列
    merged_df = merged_df.drop(columns='ID')

    # 保存合并后的数据到新的 CSV 文件
    output_file = f'{output_folder}/ped_correction_new.csv'
    merged_df.to_csv(output_file, index=False)
    print(f'新文件已生成: {output_file}')
    print(f"ped_correction_new.csv文件中共有{len(merged_df)}行数据（不算列名）")
    
    # 修改 ped_error.csv 文件，只保留表头
    if os.path.exists(ped_error_file):
        # 读取文件获取表头
        df_error = pd.read_csv(ped_error_file, nrows=0, encoding='gbk')
        # 创建只有表头的空DataFrame
        empty_df = pd.DataFrame(columns=df_error.columns)
        # 保存只有表头的文件
        empty_df.to_csv(ped_error_file, index=False)
        print(f'已清空 {ped_error_file} 的内容，只保留表头')
    else:
        # 如果文件不存在，则创建只有指定表头的空文件
        empty_df = pd.DataFrame(columns=["ID", "Fa_ID", "Mo_ID", "Fa_infor", "Mo_infor", "ori_Fa_ID", "ori_Mo_ID"])
        empty_df.to_csv(ped_error_file, index=False)
        print(f'已创建新的 {ped_error_file} 文件，只包含指定表头')

except FileNotFoundError:
    print('错误: 未找到指定的文件，请检查文件路径是否正确。')
except Exception as e:
    print(f'错误: 发生了一个未知错误: {e}')