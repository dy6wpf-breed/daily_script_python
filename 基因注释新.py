from Bio import Entrez
import pandas as pd
from xml.etree import ElementTree as ET
from io import StringIO

Entrez.email = "dy6wpf@gmail.com"

def fetch_gene_generif(gene_id):
    try:
        # 获取XML格式的基因记录
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="docsum", retmode="xml")
        xml_data = handle.read()
        handle.close()

        # 将字节序列解码为字符串
        xml_string = xml_data.decode('utf-8')

        # 使用StringIO来模拟文件操作
        xml_file = StringIO(xml_string)

        # 解析XML数据
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # 查找GeneRIF字段
        gene_rif_texts = []
        for gene_rif in root.findall('.//GeneRIF'):
            # 确保Description标签存在并提取文本内容
            description = gene_rif.find('Description')
            if description is not None:
                gene_rif_texts.append(description.text)

        # 返回GeneRIF文本，如果没有找到，则返回默认信息
        return '\n'.join(gene_rif_texts) if gene_rif_texts else "No GeneRIF available"
    except Exception as e:
        return str(e)

# 加载基因ID列表
file_path = r'E:\desk\基因注释1.xlsx'
df = pd.read_excel(file_path)

# 从Excel文件中获取基因ID
gene_ids = df['NCBI GeneID'].tolist()

# 获取每个基因的GeneRIF信息
gene_generif_list = []
for gene_id in gene_ids:
    gene_rif = fetch_gene_generif(gene_id)
    gene_generif_list.append(gene_rif)

# 将GeneRIF信息添加到DataFrame
df['GeneRIF'] = gene_generif_list

# 保存更新后的DataFrame到新的Excel文件
output_file_path = r'E:\desk\基因注释1_更新_GeneRIF.xlsx'
df.to_excel(output_file_path, index=False)

print(f"Updated file saved to: {output_file_path}")