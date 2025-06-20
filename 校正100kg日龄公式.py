import pandas as pd


def calculate_corrected_values(df):
    # 筛选大白猪公猪数据
    df = df[(df['品种品系'] == '大白') & (df['性别'] == '公')]

    # 计算达100Kg体重日龄
    A = 50.775
    df['校正100公斤体重日龄'] = df.apply(
        lambda row: row['结测日龄']+(100 - row['结测体重'])*(row['结测日龄'] - A)/row['结测体重'], axis=1)

    # 计算100Kg体重活体背膘厚
    B = -7.277
    df['校正100公斤背膘厚'] = df.apply(
        lambda row: row['结测背膘厚']+(100 - row['结测体重'])*row['结测背膘厚']/(row['结测体重'] - B), axis=1)

    # 提取需要的字段
    result = df[['个体号', '校正100公斤体重日龄', '校正100公斤背膘厚']]

    return result


if __name__ == "__main__":
    # 读取Excel文件，这里假设文件名为供精公猪系谱.XLSX，需根据实际情况修改
    try:
        excel_file =r"F:\旧电脑文件\王鹏飞\WeChat Files\WeChat Files\dy6wpf\FileStorage\File\2025-02\供精公猪系谱.XLSX"
        df = pd.read_excel(excel_file)
        result = calculate_corrected_values(df)
        print(result)
    except FileNotFoundError:
        print(f"找不到文件：{excel_file}，请检查文件名和路径是否正确。")
    except Exception as e:
        print(f"发生错误：{e}")
