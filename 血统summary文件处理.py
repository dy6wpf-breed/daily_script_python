import pandas as pd

# 读取 breed_summary.xlsx 文件
try:
    breed_summary = pd.read_excel("breed_summary.xlsx")
except FileNotFoundError:
    print("错误：找不到 breed_summary.xlsx 文件")
    exit()

# 创建样本编号到类型的字典
breed_dict = {}

# 遍历每一行数据
for index, row in breed_summary.iterrows():
    sample_ids = row["样本编号"].split(",")
    breed_type = row["类型"]
    for sample_id in sample_ids:
        breed_dict[sample_id] = breed_type

# 读取 admixture.3.csv 文件
try:
    admixture = pd.read_csv("admixture.3.csv", header=0)
except FileNotFoundError:
    print("错误：找不到 admixture.3.csv 文件")
    exit()

# 将 admixture.3.csv 样本编号转换为字符串
admixture[admixture.columns[0]] = admixture[admixture.columns[0]].astype(str)

# 创建 breed 列
admixture["breed"] = ""

# 匹配 breed 信息
matched_count = 0
breed_counts = {}
for index, row in admixture.iterrows():
    sample_id = str(row[admixture.columns[0]])
    if sample_id in breed_dict:
        breed = breed_dict[sample_id]
        admixture.loc[index, "breed"] = breed
        matched_count += 1
        if breed in breed_counts:
            breed_counts[breed] += 1
        else:
            breed_counts[breed] = 1

# 保存修改后的 admixture.3.csv 文件
admixture.to_csv("admixture_processed.csv", index=False)

print("处理完成！")
print(f"总共处理了 {len(admixture)} 个样本，成功匹配了 {matched_count} 个类型。")
for breed, count in breed_counts.items():
    print(f"匹配到类型 '{breed}' 的样本数量为 {count}。")