import pandas as pd

def get_actual_type_str(series):
    """获取Series中第一个非空值的Python类型名称"""
    series_cleaned = series.dropna()
    if series_cleaned.empty:
        return "该列无非空值"
    try:
        return type(series_cleaned.iloc[0]).__name__
    except IndexError:
        return "该列无非空值"

def compare_excel_column_formats(file1_path, file2_path):
    """比较两个Excel文件中对应列的数据格式"""
    print(f"开始比较文件: '{file1_path}' 和 '{file2_path}'\n")

    try:
        df1 = pd.read_excel(file1_path)
        print(f"成功读取文件: {file1_path}")
    except FileNotFoundError:
        print(f"错误: 文件 '{file1_path}' 未找到。")
        return
    except Exception as e:
        print(f"读取文件 '{file1_path}' 时发生错误: {e}")
        return

    try:
        df2 = pd.read_excel(file2_path)
        print(f"成功读取文件: {file2_path}")
    except FileNotFoundError:
        print(f"错误: 文件 '{file2_path}' 未找到。")
        return
    except Exception as e:
        print(f"读取文件 '{file2_path}' 时发生错误: {e}")
        return

    print("\n--- 列名和顺序比较 ---")
    cols1 = df1.columns.tolist()
    cols2 = df2.columns.tolist()

    if cols1 == cols2:
        print("两个文件的列名和顺序完全一致。")
        common_cols = cols1
    else:
        print("警告: 两个文件的列名或顺序不完全一致。")
        print(f"  文件1的列: {cols1}")
        print(f"  文件2的列: {cols2}")
        common_cols = [col for col in cols1 if col in cols2]
        if not common_cols:
            print("两个文件没有共同的列名，无法进行比较。")
            return
        print(f"将仅比较以下共同的列: {common_cols}")
        # 为了按文件1的顺序比较共同列
        common_cols = [col for col in cols1 if col in common_cols]


    print("\n--- 各列数据格式比较 ---")
    for col_name in common_cols:
        print(f"\n比较列: '{col_name}'")
        
        dtype1 = df1[col_name].dtype
        dtype2 = df2[col_name].dtype
        
        print(f"  文件1 ('{os.path.basename(file1_path)}'):")
        print(f"    Pandas dtype: {dtype1}")
        actual_type1 = get_actual_type_str(df1[col_name])
        print(f"    首个非空值Python类型: {actual_type1}")

        print(f"  文件2 ('{os.path.basename(file2_path)}'):")
        print(f"    Pandas dtype: {dtype2}")
        actual_type2 = get_actual_type_str(df2[col_name])
        print(f"    首个非空值Python类型: {actual_type2}")

        if dtype1 != dtype2:
            print(f"  !! 差异警告: Pandas dtype 不同 ({dtype1} vs {dtype2})")
        elif actual_type1 != actual_type2 and not (actual_type1 == "该列无非空值" and actual_type2 == "该列无非空值"):
             # 如果dtype相同，但实际类型不同（且都不是因为空列导致）
            print(f"  !! 差异警告: 首个非空值的实际Python类型不同 ({actual_type1} vs {actual_type2})，尽管Pandas dtype相同。")
        else:
            print(f"  数据格式一致。")
            
    print("\n比较完成。")

if __name__ == "__main__":
    import os
    # 请确保文件路径正确
    # 假设脚本和Excel文件都在 'y:\种质资源数据管理平台2024\9-自贡德康\' 目录下
    base_dir = r"y:\种质资源数据管理平台2024\9-自贡德康"
    file1 = os.path.join(base_dir, "1-样本基因型对应关系表1.xlsx")
    file2 = os.path.join(base_dir, "1-样本基因型对应关系表2.xlsx")
    
    compare_excel_column_formats(file1, file2)