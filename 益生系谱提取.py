import pandas as pd
import os

# 指定 Excel 文件的路径
# file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【2】黑龙江益生-基因型数据TOP链要先转至GGP50K\【100】 KPSNY132025003C_KPS202505088_黑龙江益生种猪繁育有限公司2543_唐友_40个-20250526\KPSNY132025003C_唐友_黑龙江益生种猪繁育有限公司2543_40个_中芯一号PLUS_report_20250526\1.Sampleinfo\2543批次系谱.xlsx"

# 提示用户手动输入文件路径
file_path_input = input("请输入 Excel 文件路径: ")

# 移除路径两端的双引号
file_path = file_path_input.strip('"')

try:
    # 读取 Excel 文件
    df = pd.read_excel(file_path)
    # 去除列名中的空格
    df.columns = df.columns.str.strip()

    # 提取所需的字段
    columns_to_extract = ["个体号", "品种品系", "出生日期", "性别", "父亲", "母亲"]
    # 使用 .copy() 确保 new_df 是一个副本
    new_df = df[columns_to_extract].copy()

    # 根据出生日期的年份生成 year 字段
    new_df["year"] = pd.to_datetime(new_df["出生日期"]).dt.year

    # 将性别中的"公"改成 M，"母"改成 F
    new_df["性别"] = new_df["性别"].replace({"公": "M", "母": "F"})

    # 调整列的顺序
    new_df = new_df[["个体号", "品种品系", "出生日期", "year", "性别", "父亲", "母亲"]]

    # 获取文件所在的文件夹路径
    folder_path = os.path.dirname(file_path)
    # 生成新表格的路径
    new_file_path = os.path.join(folder_path, "新表格.xlsx")

    # 将新表格保存为 Excel 文件
    new_df.to_excel(new_file_path, index=False)
    print(f"新表格已保存到 {new_file_path}")
except FileNotFoundError:
    print(f"未找到文件: {file_path}")
except KeyError as e:
    print(f"提取列时出错，指定的列名可能不存在: {e}")
except Exception as e:
    print(f"处理文件时出现未知错误: {e}")
    