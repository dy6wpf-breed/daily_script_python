import pandas as pd
import os
import glob
from datetime import datetime

def find_chip_position_files(root_folder):
    """
    遍历文件夹及其子文件夹，找到文件名包含"芯片位置对应表"的Excel文件
    
    参数:
        root_folder: 根文件夹路径
        
    返回:
        包含所有匹配文件路径的列表
    """
    print(f"正在搜索文件夹: {root_folder}")
    matched_files = []
    
    # 使用glob递归搜索所有Excel文件
    excel_patterns = ['*.xlsx', '*.xls']
    for pattern in excel_patterns:
        search_pattern = os.path.join(root_folder, '**', pattern)
        for file_path in glob.glob(search_pattern, recursive=True):
            # 检查文件名是否包含"芯片位置对应表"
            if "芯片位置对应表" in os.path.basename(file_path):
                matched_files.append(file_path)
                print(f"找到匹配文件: {file_path}")
    
    return matched_files

def extract_data_from_file(file_path):
    """
    从Excel文件中提取指定列的数据
    
    参数:
        file_path: Excel文件路径
        
    返回:
        包含提取数据的DataFrame
    """
    print(f"正在处理文件: {file_path}")
    
    try:
        # 读取Excel文件，从第5行开始（第5行是列名）
        df = pd.read_excel(file_path, header=4)
        
        # 获取文件创建时间
        file_creation_time = datetime.fromtimestamp(os.path.getctime(file_path))
        file_creation_date = file_creation_time.strftime('%Y/%m/%d')
        file_creation_year = file_creation_time.strftime('%Y')
        
        # 需要提取的列
        columns_to_extract = [
            "任务单号", "合同号", "客户单位", "Sample_ID", 
            "实验室编号", "客户姓名", "本次上机数目", "Sample_Name"
        ]
        
        # 检查列是否存在
        available_columns = []
        for col in columns_to_extract:
            if col in df.columns:
                available_columns.append(col)
            else:
                print(f"警告: 列 '{col}' 在文件 {file_path} 中不存在")
        
        # 如果没有找到任何需要的列，返回空DataFrame
        if not available_columns:
            print(f"错误: 文件 {file_path} 中没有找到任何需要提取的列")
            return pd.DataFrame()
        
        # 提取可用的列
        extracted_df = df[available_columns].copy()
        
        # 添加文件创建日期和年份
        extracted_df["下机日期"] = file_creation_date
        extracted_df["检测年份"] = file_creation_year
        
        # 添加文件路径信息用于调试
        extracted_df["源文件"] = os.path.basename(file_path)
        
        return extracted_df
    
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {str(e)}")
        return pd.DataFrame()

def map_columns_to_new_format(df):
    """
    将提取的数据映射到新的表格格式
    
    参数:
        df: 包含提取数据的DataFrame
        
    返回:
        映射后的DataFrame
    """
    # 定义列映射关系
    column_mapping = {
        "任务单号": "任务单号",
        "合同号": "合同编号",
        "客户单位": "客户单位",
        "Sample_ID": "Sample_ID",
        "实验室编号": "实验室编号",
        "客户姓名": "客户姓名",
        "下机日期": "下机日期",
        "检测年份": "检测年份",
        "本次上机数目": "样本数",
        "Sample_Name": "个体ID"
    }
    
    # 创建新的DataFrame，包含所有目标列
    new_columns = [
        "任务单号", "合同编号", "客户单位", "品种", "品系", "数据类型",
        "Sample_ID", "实验室编号", "客户姓名", "下机日期", "检测年份",
        "批次", "样本数", "个体ID"
    ]
    result_df = pd.DataFrame(columns=new_columns)
    
    # 复制数据到新的DataFrame
    for old_col, new_col in column_mapping.items():
        if old_col in df.columns and new_col in result_df.columns:
            result_df[new_col] = df[old_col]
    
    return result_df

def remove_duplicates(df):
    """
    根据Sample_ID列删除重复行
    
    参数:
        df: 包含数据的DataFrame
        
    返回:
        去重后的DataFrame
    """
    if "Sample_ID" not in df.columns:
        print("警告: 数据中没有Sample_ID列，无法进行去重")
        return df
    
    # 记录去重前的行数
    before_count = len(df)
    
    # 检查是否有重复的Sample_ID
    duplicate_ids = df[df["Sample_ID"].duplicated()]["Sample_ID"].unique()
    if len(duplicate_ids) > 0:
        print(f"发现 {len(duplicate_ids)} 个重复的Sample_ID: {', '.join(map(str, duplicate_ids[:5]))}")
        if len(duplicate_ids) > 5:
            print(f"... 以及其他 {len(duplicate_ids) - 5} 个")
    
    # 去除重复行，保留第一次出现的行
    df_no_duplicates = df.drop_duplicates(subset=["Sample_ID"], keep="first")
    
    # 记录去重后的行数
    after_count = len(df_no_duplicates)
    removed_count = before_count - after_count
    
    print(f"去重前行数: {before_count}")
    print(f"去重后行数: {after_count}")
    print(f"已删除 {removed_count} 行重复数据")
    
    return df_no_duplicates

def main():
    # 获取当前脚本所在目录作为根目录
    root_folder = os.path.dirname(os.path.abspath(__file__))
    
    # 询问用户是否要使用其他目录
    use_current_dir = input(f"是否使用当前目录 ({root_folder}) 作为搜索根目录? (y/n): ").strip().lower()
    
    if use_current_dir != 'y':
        root_folder = input("请输入要搜索的根目录路径: ").strip()
        # 检查输入的路径是否存在
        if not os.path.exists(root_folder):
            print(f"错误: 目录 '{root_folder}' 不存在!")
            return
    
    # 查找匹配的文件
    matched_files = find_chip_position_files(root_folder)
    
    if not matched_files:
        print("未找到任何包含'芯片位置对应表'的Excel文件")
        return
    
    print(f"共找到 {len(matched_files)} 个匹配的文件")
    
    # 处理每个文件并合并数据
    all_data = []
    for file_path in matched_files:
        extracted_data = extract_data_from_file(file_path)
        if not extracted_data.empty:
            all_data.append(extracted_data)
    
    if not all_data:
        print("未能从任何文件中提取有效数据")
        return
    
    # 合并所有提取的数据
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # 映射到新的表格格式
    result_df = map_columns_to_new_format(combined_df)
    
    # 对Sample_ID列进行去重
    result_df = remove_duplicates(result_df)
    
    # 保存结果
    output_file = os.path.join(root_folder, "芯片位置对应表汇总.xlsx")
    result_df.to_excel(output_file, index=False)
    
    print(f"数据提取、去重和汇总完成! 已保存到: {output_file}")
    print(f"共处理了 {len(matched_files)} 个文件，提取了 {len(result_df)} 条记录")

if __name__ == "__main__":
    main()