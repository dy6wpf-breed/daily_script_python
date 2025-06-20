import pandas as pd
from Bio import Entrez
import time

def get_gene_id(gene_name, organism="Bos taurus"):
    Entrez.email = "bjlmwpf@163.com"  # 请替换为您的邮箱
    handle = Entrez.esearch(db="gene", term=f"{gene_name}[Gene Name] AND {organism}[Organism]")
    record = Entrez.read(handle)

    if record["Count"] == "0":
        return f"没有找到 {gene_name} 的基因ID"
    else:
        return record["IdList"][0]


# 读取Excel文件
file_path = r"E:\desk\基因名称.xlsx"
df = pd.read_excel(file_path)

# 获取第一列的所有基因名称
gene_names = df.iloc[:, 0].tolist()

# 创建一个字典来存储结果
results = {}

# 查询每个基因名称的Gene ID
for gene_name in gene_names:
    if pd.notna(gene_name):  # 检查是否为空值
        gene_id = get_gene_id(gene_name)
        results[gene_name] = gene_id
        print(f"基因名称: {gene_name}, Gene ID: {gene_id}")
        time.sleep(1)  # 添加延迟以避免过快请求NCBI服务器

# 将结果保存到新的Excel文件
result_df = pd.DataFrame(list(results.items()), columns=['基因名称', 'Gene ID'])
result_df.to_excel(r"E:\desk\基因ID结果.xlsx", index=False)

print("查询完成，结果已保存到 'E:\desk\基因ID结果.xlsx'")