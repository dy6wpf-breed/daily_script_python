import pandas as pd
import os
import re
from datetime import datetime

def reorganize_breeding_data(source_file, template_file, output_file=None):
    """
    根据模板文件的结构整理繁殖表型数据
    
    参数:
        source_file: 源数据文件路径 (繁殖表型.xlsx)
        template_file: 模板文件路径 (多次动物表型数据表.xlsx)
        output_file: 输出文件路径，默认为"整理后的繁殖表型数据.xlsx"
    """
    if output_file is None:
        output_file = "整理后的繁殖表型数据.xlsx"
    
    print(f"正在读取源数据文件: {source_file}")
    source_data = pd.read_excel(source_file)
    
    print(f"正在读取模板文件: {template_file}")
    template_data = pd.read_excel(template_file)
    
    # 获取模板的列名和结构
    template_columns = template_data.columns.tolist()
    
    print(f"模板文件包含以下列: {template_columns}")
    print(f"源数据文件包含以下列: {source_data.columns.tolist()}")
    
    # 根据提供的对应关系进行映射
    column_mapping = {
        "客户单位": "黑龙江益生",  # 固定值
        "场": "本地猪出生猪场",    # 固定值
        "个体ID": "母猪号",
        "管理号": "场内编号",
        "性别": "性别",
        "配种日期": None,
        "胎次": "胎次",
        "与配公猪": "与配公猪",
        "分娩日期": "分娩日期",
        "繁殖性状-总仔": "总仔",
        "繁殖性状-活仔": "活仔",
        "繁殖性状-健仔": "健仔",
        "断奶日期": "断奶日期",
        "繁殖性状-断奶头数": "断奶头数",
        "繁殖性状-哺乳头数": None,
        "21日称重日期": "21日称重日期",
        "繁殖性状-21天窝重": "21天窝重",
        "繁殖性状-21日称重头数": "21日称重头数",
        "繁殖性状-哺乳期成活率": None,
        "繁殖性状-校正21日窝重": "校正21日窝重",
        "断奶至首次配种天数": None,
        "妊娠天数": None,
        "哺乳天数": None
    }
    
    # 初始化结果DataFrame，确保有足够的行数
    row_count = len(source_data)
    result_data = pd.DataFrame(index=range(row_count), columns=template_columns)
    
    # 处理每一列
    for template_col in template_columns:
        # 1. 检查是否有特定映射
        if template_col in column_mapping:
            mapping_value = column_mapping[template_col]
            
            # 处理固定值的情况（客户单位和场）
            if template_col in ["客户单位", "场"] and isinstance(mapping_value, str):
                result_data[template_col] = mapping_value
            # 处理从源数据映射的列
            elif mapping_value is not None and mapping_value in source_data.columns:
                result_data[template_col] = source_data[mapping_value].values
            # 其他映射为None的情况
            else:
                result_data[template_col] = None
        # 2. 检查是否有相同列名（没有特别提到的情况）
        elif template_col in source_data.columns:
            # 如果模板列名在源数据中存在，直接取值
            result_data[template_col] = source_data[template_col].values
        # 3. 其他情况创建空列
        else:
            # 对于没有映射的列，创建空列
            result_data[template_col] = None
    
    # 处理日期格式（确保所有日期列格式为年/月/日）
    date_columns = ["配种日期", "分娩日期", "断奶日期", "21日称重日期"]
    for col in date_columns:
        if col in result_data.columns:
            result_data[col] = result_data[col].apply(lambda x: format_date(x) if pd.notna(x) else None)
    
    # 保存结果
    print(f"正在保存整理后的数据到: {output_file}")
    result_data.to_excel(output_file, index=False)
    print(f"数据整理完成! 已保存到: {output_file}")
    
    return result_data

def format_date(date_val):
    """格式化日期为年/月/日格式"""
    if pd.isna(date_val):
        return None
    try:
        if isinstance(date_val, str):
            # 尝试解析字符串日期
            date_obj = pd.to_datetime(date_val)
        else:
            # 已经是日期对象
            date_obj = pd.to_datetime(date_val)
        # 返回格式化的日期字符串 (年/月/日)
        return date_obj.strftime('%Y/%m/%d')
    except:
        return date_val

def main():
    # 文件路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    source_file = os.path.join(current_dir, "繁殖表型.xlsx")
    template_file = os.path.join(current_dir, "多次动物表型数据表.xlsx")
    output_file = os.path.join(current_dir, "整理后的繁殖表型数据.xlsx")
    
    # 检查文件是否存在
    if not os.path.exists(source_file):
        print(f"错误: 源数据文件 '{source_file}' 不存在!")
        return
    
    if not os.path.exists(template_file):
        print(f"错误: 模板文件 '{template_file}' 不存在!")
        return
    
    # 执行数据整理
    reorganize_breeding_data(source_file, template_file, output_file)

if __name__ == "__main__":
    main()