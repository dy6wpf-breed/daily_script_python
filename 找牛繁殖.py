import pandas as pd
import re
from tqdm import tqdm
import requests
import time

def parse_gff(filename):
    print("正在读取GFF文件...")
    gene_info = []
    
    reproduction_keywords = [
        'reproduction', 'fertility', 'ovarian', 'testis', 'sperm', 'oocyte',
        'estrus', 'pregnancy', 'embryo', 'placenta', 'uterus', 'ovulation',
        'spermatogenesis', 'oogenesis', 'follicle', 'hormone', 'gonad'
    ]

    keyword_pattern = re.compile('|'.join(reproduction_keywords), re.IGNORECASE)

    try:
        with open(filename, 'r') as f:
            for line in tqdm(f, desc="处理GFF文件"):
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) != 9:
                    continue
                if parts[2] == 'gene':
                    attributes = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
                    description = attributes.get('description', '')
                    if keyword_pattern.search(description):
                        gene_info.append({
                            'Gene Name': attributes.get('Name', ''),
                            'Gene ID': attributes.get('ID', ''),
                            'Chromosome': parts[0],
                            'Start': parts[3],
                            'End': parts[4],
                            'Strand': parts[6],
                            'Description': description
                        })
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        return []

    return gene_info

def get_gene_description(gene_id):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch_url = f"{base_url}esearch.fcgi?db=gene&term={gene_id}[Gene Name] AND bos taurus[Organism]&retmode=json"
    
    try:
        response = requests.get(esearch_url)
        data = response.json()
        if 'esearchresult' in data and 'idlist' in data['esearchresult'] and data['esearchresult']['idlist']:
            ncbi_id = data['esearchresult']['idlist'][0]
            esummary_url = f"{base_url}esummary.fcgi?db=gene&id={ncbi_id}&retmode=json"
            summary_response = requests.get(esummary_url)
            summary_data = summary_response.json()
            if 'result' in summary_data and ncbi_id in summary_data['result']:
                return summary_data['result'][ncbi_id]['description']
    except Exception as e:
        print(f"获取基因 {gene_id} 的描述时发生错误: {e}")
    
    return "无法获取描述"

def main():
    # 使用函数解析GFF文件
    gene_info = parse_gff(r"E:\desk\Bos_taurus.ARS-UCD1.3.112.gff3")

    if gene_info:
        # 创建DataFrame
        df = pd.DataFrame(gene_info)

        print("正在获取基因功能描述...")
        # 为每个基因获取功能描述
        for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="获取基因描述"):
            gene_id = row['Gene Name']
            df.at[index, 'Function Description'] = get_gene_description(gene_id)
            time.sleep(0.5)  # 为了避免对NCBI服务器发送过多请求，我们在每次请求之间暂停0.5秒

        print("正在导出数据到Excel...")
        # 导出到Excel
        df.to_excel('bovine_reproduction_genes_with_function.xlsx', index=False)

        print(f"找到 {len(gene_info)} 个可能与繁殖相关的基因。数据已成功导出到 bovine_reproduction_genes_with_function.xlsx")
    else:
        print("未找到相关基因信息。")

if __name__ == "__main__":
    main()