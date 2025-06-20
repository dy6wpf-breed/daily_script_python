import pandas as pd
import os

# 提示用户输入文件路径
file_path = input("请输入 Excel 文件的路径: ")
# 移除路径两端的双引号
file_path = file_path.strip('"')

try:
    # 读取 Excel 文件
    df = pd.read_excel(file_path)
    # 去除列名中的空格
    df.columns = df.columns.str.strip()

    # 提取所需的字段，将 "SIRE" 替换为 "父亲ID"，"DAM" 替换为 "母亲ID"
    columns_to_extract = ["种猪ID", "管理号", "品种", "品系", "性别", "出生日期", "父亲ID", "母亲ID"]
    # 使用 .copy() 确保 new_df 是一个副本
    new_df = df[columns_to_extract].copy()

    # 重命名列名，将 "父亲ID" 改为 "SIRE"，"母亲ID" 改为 "DAM"
    new_df.rename(columns={"父亲ID": "SIRE", "母亲ID": "DAM"}, inplace=True)

    # 根据出生日期的年份生成 YEAR 字段
    new_df["YEAR"] = pd.to_datetime(new_df["出生日期"]).dt.year

    # 将性别中的"公"改成 M，"母"改成 F
    new_df["性别"] = new_df["性别"].replace({"公": "M", "母": "F"})

    # 调整列的顺序
    new_df = new_df[["种猪ID", "管理号", "品种", "品系", "性别", "出生日期", "YEAR", "SIRE", "DAM"]]

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