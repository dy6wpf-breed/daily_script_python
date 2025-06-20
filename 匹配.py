import pandas as pd
import gzip
import os
from deep_translator import GoogleTranslator
import time
from Bio import Entrez

# 确认文件路径
file_path = r'E:\desk\基因注释1.xlsx'
gz_file_path = r'E:\desk\generifs_basic.gz'  # 修改为文件的实际路径

# 检查文件是否存在
if not os.path.exists(file_path):
    raise FileNotFoundError(f"Excel文件未找到：{file_path}")
if not os.path.exists(gz_file_path):
    raise FileNotFoundError(f"gz文件未找到：{gz_file_path}")

# 读取Excel文件
df = pd.read_excel(file_path)

# 读取并解析generifs_basic.gz文件，将其转换为DataFrame
def read_gz_file_to_df(gz_file_path):
    data = []
    with gzip.open(gz_file_path, 'rt') as f:
        # 逐行读取文件内容并解析
        for line in f:
            if line.startswith('#'):
                continue  # 跳过标题行
            fields = line.strip().split('\t')
            if len(fields) == 5:
                data.append(fields)
    # 将解析的数据转换为DataFrame
    columns = ['TaxID', 'GeneID', 'PMID_list', 'timestamp', 'GeneRIF']
    return pd.DataFrame(data, columns=columns)

# 从generifs_basic.gz文件中提取GeneRIF信息
gene_rifs_df = read_gz_file_to_df(gz_file_path)

# 打印DataFrame的前几行以检查结构
print("gene_rifs_df前几行:\n", gene_rifs_df.head())

# 确认DataFrame中的列名
print("gene_rifs_df列名:\n", gene_rifs_df.columns)

# 初始化翻译器
translator = GoogleTranslator(source='en', target='zh-CN')

# 设置Entrez邮箱
Entrez.email = "bjlmwpf@163.com"

# 翻译GeneRIF信息
def translate_rif(text, retries=3):
    if text == "未找到GeneRIF":
        return text
    while retries > 0:
        try:
            translation = translator.translate(text)
            return translation
        except Exception as e:
            retries -= 1
            time.sleep(1)  # 等待1秒后重试
    return text  # 返回原始英文内容

# 从NCBI获取基因注释
def fetch_gene_annotation(gene_id):
    try:
        handle = Entrez.esummary(db="gene", id=gene_id)
        record = Entrez.read(handle)
        handle.close()
        if record and 'Summary' in record[0]:
            return record[0]['Summary']
    except Exception as e:
        print(f"Error fetching annotation for GeneID {gene_id}: {e}")
    return "未找到GeneRIF"

# 合并GeneRIF信息并翻译
def get_translated_gene_rif(gene_id):
    rifs = gene_rifs_df[gene_rifs_df['GeneID'] == str(gene_id)]['GeneRIF'].tolist()
    if rifs:
        concatenated_rifs = " ".join(rifs)
        translated_rif = translate_rif(concatenated_rifs)
        return translated_rif
    else:
        annotation = fetch_gene_annotation(gene_id)
        return translate_rif(annotation)

# 创建一个新的DataFrame来存储翻译结果
translated_results = []

for gene_id in df['NCBI GeneID']:  # 修改为实际的列名
    translated_rif = get_translated_gene_rif(gene_id)
    translated_results.append({
        'GeneID': gene_id,
        'TranslatedGeneRIF': translated_rif
    })

translated_df = pd.DataFrame(translated_results)

# 输出结果到新的Excel文件
output_file_path = r'E:\desk\翻译结果.xlsx'
translated_df.to_excel(output_file_path, index=False)

print(f"翻译结果已保存到：{output_file_path}")