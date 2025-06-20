import pandas as pd
from Bio import Entrez

# 设置NCBI的电子邮件地址
Entrez.email = "dy6wpf@gmail.com"

def get_gene_location(gene_id, species="goat"):
    try:
        search_query = f"{gene_id}[Gene Name] AND {species}[Orgn]"
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=10)
        search_results = Entrez.read(handle)
        handle.close()

        if not search_results['IdList']:
            print(f"No results found for gene ID: {gene_id}")
            return

        for id in search_results['IdList']:
            handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="binary")
            record = Entrez.read(handle)
            handle.close()

            for feature in record['GBSeq_feature-table']:
                if feature['FTName'] == 'gene':
                    gene_location = feature['FTLocation']
                    print(f"Gene ID: {gene_id}")
                    print(f"Species: {species}")
                    print(f"Chromosome: {record['GBSeq_source']}")
                    print(f"Location: {gene_location}")
                    return
        print(f"Gene location not found for ID: {id}")
    except Exception as e:
        print(f"An error occurred: {e}")

def read_gene_ids_from_excel(file_path):
    try:
        df = pd.read_excel(file_path)
        gene_ids = df['NCBI GeneID'].drop_duplicates().tolist()
        return gene_ids
    except Exception as e:
        print(f"Failed to read Excel file: {e}")
        return []

# 示例使用
file_path = r"E:\desk\基因注释1.xlsx"  # 替换为你的Excel文件路径
gene_ids = read_gene_ids_from_excel(file_path)
for gene_id in gene_ids:
    get_gene_location(str(gene_id))