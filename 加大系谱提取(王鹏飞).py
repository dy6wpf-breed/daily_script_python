import pandas as pd
import os

# 指定 Excel 文件的路径
# file_path = r"Z:\猪中芯一号V1PLUS\【1】GS分析\【5】江西加大\【14】KPSNY032024027C_KPS202505108_江西加大集团有限公司_帅亮_96个中芯一号mini-20250526\加大系谱.xlsx"

# 手动输入文件路径并移除双引号
file_path = input("请输入 Excel 文件路径: ").strip('"')

try:
    # 读取 Excel 文件
    df = pd.read_excel(file_path)
    # 去除列名中的空格
    df.columns = df.columns.str.strip()

    # 提取所需的字段
    columns_to_extract = ["个体号", "耳缺号（耳牌号）", "场内编号", "本地猪出生猪场", "品种品系", "出生日期", "性别", "父亲", "母亲"]
    # 使用 .copy() 确保 new_df 是一个副本
    new_df = df[columns_to_extract].copy()

    # 根据出生日期的年份生成 YEAR 字段
    new_df["YEAR"] = pd.to_datetime(new_df["出生日期"]).dt.year

    # 将性别中的"公"改成 M，"母"改成 F
    new_df["性别"] = new_df["性别"].replace({"公": "M", "母": "F"})

    # 调整列的顺序
    new_df = new_df[["个体号", "耳缺号（耳牌号）", "场内编号", "本地猪出生猪场", "品种品系", "出生日期", "YEAR", "性别", "父亲", "母亲"]]

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