import pandas as pd
import os

print(f"pandas version: {pd.__version__}")  # 添加pandas版本信息

def calculate_and_standardize_ebv(file_path, output_file_path):
    """
    计算综合育种值并进行标准化。

    Args:
        file_path (str): Excel文件路径。
        output_file_path (str): 输出Excel文件路径。
    """
    try:
        print("开始执行脚本...")  # 添加开始信息

        # 读取Excel文件
        print("正在读取Excel文件...")  # 添加读取信息
        df_all = pd.read_excel(file_path, sheet_name="全部分析个体")
        df_current = pd.read_excel(file_path, sheet_name="本次分析个体")

        # 假设列名，如果列名不正确，会抛出异常
        总仔数_blup_col = "健仔数_blup"
        矫正背膘厚_blup_col = "矫正背膘厚_blup"
        矫正百公斤日龄_blup_col = "矫正百公斤日龄_blup"

        # 标准化每个性状
        print("正在标准化每个性状...")  # 添加标准化信息

        # 全部分析个体
        mean_tnb_all = df_all[总仔数_blup_col].mean()
        std_tnb_all = df_all[总仔数_blup_col].std()
        # df_all["标准化健仔数"] = 100 + 25 * (df_all[总仔数_blup_col] - mean_tnb_all) / std_tnb_all # Old single-trait standardization

        mean_bf_all = df_all[矫正背膘厚_blup_col].mean()
        std_bf_all = df_all[矫正背膘厚_blup_col].std()
        # df_all["标准化矫正背膘厚_blup"] = 100 + 25 * (df_all[矫正背膘厚_blup_col] - mean_bf_all) / std_bf_all # Old single-trait standardization

        mean_day_all = df_all[矫正百公斤日龄_blup_col].mean()
        std_day_all = df_all[矫正百公斤日龄_blup_col].std()
        # df_all["标准化矫正百公斤日龄_blup"] = 100 + 25 * (df_all[矫正百公斤日龄_blup_col] - mean_day_all) / std_day_all # Old single-trait standardization

        # 本次分析个体
        mean_tnb_current = df_current[总仔数_blup_col].mean()
        std_tnb_current = df_current[总仔数_blup_col].std()
        # df_current["标准化健仔数"] = 100 + 25 * (df_current[总仔数_blup_col] - mean_tnb_current) / std_tnb_current # Old single-trait standardization

        mean_bf_current = df_current[矫正背膘厚_blup_col].mean()
        std_bf_current = df_current[矫正背膘厚_blup_col].std()
        # df_current["标准化矫正背膘厚_blup"] = 100 + 25 * (df_current[矫正背膘厚_blup_col] - mean_bf_current) / std_bf_current # Old single-trait standardization

        mean_day_current = df_current[矫正百公斤日龄_blup_col].mean()
        std_day_current = df_current[矫正百公斤日龄_blup_col].std()
        # df_current["标准化矫正百公斤日龄_blup"] = 100 + 25 * (df_current[矫正百公斤日龄_blup_col] - mean_day_current) / std_day_current # Old single-trait standardization

        # 计算综合育种值中间指数 I
        print("正在计算综合育种值中间指数...")
        df_all["综合育种值_I"] = (
            -0.25 * (df_all[矫正百公斤日龄_blup_col] / std_day_all)
            - 0.1 * (df_all[矫正背膘厚_blup_col] / std_bf_all)
            + 0.65 * (df_all[总仔数_blup_col] / std_tnb_all)
        )

        df_current["综合育种值_I"] = (
            -0.25 * (df_current[矫正百公斤日龄_blup_col] / std_day_current)
            - 0.1 * (df_current[矫正背膘厚_blup_col] / std_bf_current)
            + 0.65 * (df_current[总仔数_blup_col] / std_tnb_current)
        )

        # 对中间指数 I 进行标准化
        print("正在标准化综合育种值中间指数...")
        mean_I_all = df_all["综合育种值_I"].mean()
        std_I_all = df_all["综合育种值_I"].std()
        df_all["综合育种值"] = 100 + 25 * (df_all["综合育种值_I"] - mean_I_all) / std_I_all

        mean_I_current = df_current["综合育种值_I"].mean()
        std_I_current = df_current["综合育种值_I"].std()
        df_current["综合育种值"] = 100 + 25 * (df_current["综合育种值_I"] - mean_I_current) / std_I_current

        # 删除中间列 综合育种值_I
        df_all = df_all.drop(columns=["综合育种值_I"])
        df_current = df_current.drop(columns=["综合育种值_I"])

        # 检查文件是否存在，如果存在，则尝试删除它
        if os.path.exists(output_file_path):
            print(f"文件 {output_file_path} 存在，尝试删除...")  # 添加删除文件信息
            os.remove(output_file_path)

        # 将结果写入Excel文件
        print("正在写入Excel文件...")  # 添加写入信息
        with pd.ExcelWriter(output_file_path) as writer:
            df_all.to_excel(writer, sheet_name="全部分析个体", index=False)
            df_current.to_excel(writer, sheet_name="本次分析个体", index=False)

        print(f"综合育种值已计算并标准化，结果已保存到 {output_file_path}")  # 添加完成信息

    except KeyError as e:
        print(f"列名错误：{e}。请检查Excel文件中的列名是否正确。")
    except Exception as e:
        print(f"发生错误：{e}")

# 设置文件路径
    # 从用户输入获取文件路径并移除双引号
file_path = input("请输入Excel文件路径：").strip('\"')
# 获取输入文件所在的目录，并构建输出文件路径
output_directory = os.path.dirname(file_path)
output_file_path = os.path.join(output_directory, "BLUP_total_calculated.xlsx")

# 调用函数
calculate_and_standardize_ebv(file_path, output_file_path)