import pandas as pd
import os
import re
from datetime import datetime

def reorganize_data(source_file, template_file, output_file=None):
    """
    根据模板文件的结构整理源文件中的数据
    
    参数:
        source_file: 源数据文件路径 (表型.xlsx)
        template_file: 模板文件路径 (表型模板.xlsx)
        output_file: 输出文件路径，默认为"整理后的表型数据.xlsx"
    """
    if output_file is None:
        output_file = "整理后的表型数据.xlsx"
    
    print(f"正在读取源数据文件: {source_file}")
    source_data = pd.read_excel(source_file)
    
    print(f"正在读取模板文件: {template_file}")
    template_data = pd.read_excel(template_file)
    
    # 获取模板的列名和结构
    template_columns = template_data.columns.tolist()
    
    print(f"模板文件包含以下列: {template_columns}")
    print(f"源数据文件包含以下列: {source_data.columns.tolist()}")
    
    # 创建结果DataFrame
    result_data = pd.DataFrame(columns=template_columns)
    
    # 根据提供的对应关系进行映射
    column_mapping = {
        "客户单位": "黑龙江益生",  # 固定值
        "场": "本地猪出生猪场",    # 固定值
        "个体ID": "个体号",
        "品种": None,  # 需要从个体号提取前两位字母
        "品系": "丹系",  # 固定值
        "同窝总仔数": "同窝仔猪数",
        "左乳头数": "出生左乳头数",
        "右乳头数": "出生右乳头数",
        "达100kg体重日龄": "校正100公斤体重日龄",
        "达100kg背膘厚": "校正100公斤背膘厚(B)"
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
            if template_col in ["客户单位", "场", "品系"]:
                result_data[template_col] = mapping_value
            # 处理从源数据映射的列
            elif mapping_value is not None and mapping_value in source_data.columns:
                result_data[template_col] = source_data[mapping_value].values
            # 其他映射为None的情况
            else:
                result_data[template_col] = None
        # 2. 特殊处理品种列
        elif template_col == "品种":
            # 特殊处理：从个体号提取前两位字母
            if "个体号" in source_data.columns:
                def extract_letters(id_str):
                    if pd.isna(id_str):
                        return None
                    # 将输入转换为字符串
                    id_str = str(id_str)
                    # 使用正则表达式提取字母部分
                    letters = re.findall(r'[A-Za-z]+', id_str)
                    if letters:
                        # 返回第一组字母的前两位（如果有的话）
                        return letters[0][:2].upper() if len(letters[0]) >= 2 else letters[0].upper()
                    return None
                
                result_data[template_col] = source_data["个体号"].apply(extract_letters)
            else:
                result_data[template_col] = None
        # 3. 特殊处理出生日期列，确保格式为年/月/日
        elif template_col == "出生日期":
            if template_col in source_data.columns:
                # 处理日期格式
                def format_date(date_val):
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
                
                result_data[template_col] = source_data[template_col].apply(format_date)
            else:
                result_data[template_col] = None
        # 4. 检查是否有相同列名（没有特别提到的情况）
        elif template_col in source_data.columns:
            # 如果模板列名在源数据中存在，直接取值
            result_data[template_col] = source_data[template_col].values
        # 5. 其他情况创建空列
        else:
            # 对于没有映射的列，创建空列
            result_data[template_col] = None
    
    # 保存结果
    print(f"正在保存整理后的数据到: {output_file}")
    result_data.to_excel(output_file, index=False)
    print(f"数据整理完成! 已保存到: {output_file}")
    
    return result_data

def main():
    # 文件路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    source_file = os.path.join(current_dir, "表型.xlsx")
    template_file = os.path.join(current_dir, "表型模板.xlsx")
    output_file = os.path.join(current_dir, "整理后的表型数据.xlsx")
    
    # 检查文件是否存在
    if not os.path.exists(source_file):
        print(f"错误: 源数据文件 '{source_file}' 不存在!")
        return
    
    if not os.path.exists(template_file):
        print(f"错误: 模板文件 '{template_file}' 不存在!")
        return
    
    # 执行数据整理
    reorganize_data(source_file, template_file, output_file)

if __name__ == "__main__":
    main()