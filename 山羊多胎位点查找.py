from Bio import Entrez
import pandas as pd
from xml.etree import ElementTree as ET
from io import StringIO
import time

# 设置Entrez的电子邮件地址
Entrez.email = "dy6wpf@gmail.com"

def get_gene_id_by_name(gene_name, species):
    """通过基因名称和物种获取GeneID"""
    try:
        search_term = f"{gene_name}[Gene Name] AND {species}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_term, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        if record['IdList']:
            return record['IdList'][0]
        else:
            return None
    except Exception as e:
        return str(e)

def fetch_gene_info(gene_id):
    """根据GeneID获取基因信息"""
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="xml")
        xml_data = handle.read()
        handle.close()

        xml_string = xml_data.decode('utf-8')
        xml_file = StringIO(xml_string)
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # 获取染色体编号
        chromosome = root.findtext('.//Entrezgene_chromosome', default="Not available")
        
        # 尝试多个可能的位置来获取染色体位置信息
        chromosome_location = "Not available"
        possible_locations = [
            './/Maps/Map/Map_display-str',
            './/Entrezgene_gene/Gene-ref/Gene-ref_maploc',
            './/Entrezgene_location/Maps/Map/Map_display-str',
            './/Gene-commentary_label[text()="Chromosome"]/following-sibling::Gene-commentary_text'
        ]
        
        for location in possible_locations:
            chromosome_location = root.findtext(location)
            if chromosome_location:
                break

        # 获取基因位置
        gene_location = root.findtext('.//Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_from', default="Not available")
        gene_location += " - " + root.findtext('.//Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_to', default="Not available")

        # 获取基因注释信息
        gene_summary = root.findtext('.//Entrezgene_summary', default="Not available")

        return chromosome, chromosome_location, gene_location, gene_summary
    except Exception as e:
        return str(e), None, None, None

# 加载基因名称和物种信息
file_path = r"E:\desk\基因注释2.xlsx"
df = pd.read_excel(file_path)

# 假设Excel文件中有基因名称列和物种列
gene_names = df['基因名称'].tolist()
species_names = df['物种'].tolist()

# 获取每个基因的详细信息
gene_id_list = []
chromosome_list = []
chromosome_location_list = []
gene_location_list = []
gene_summary_list = []

for gene_name, species in zip(gene_names, species_names):
    print(f"Processing gene: {gene_name}, species: {species}")
    gene_id = get_gene_id_by_name(gene_name, species)
    if gene_id:
        gene_id_list.append(gene_id)
        chromosome, chromosome_location, gene_location, gene_summary = fetch_gene_info(gene_id)
        print(f"Gene: {gene_name}, GeneID: {gene_id}, Chromosome: {chromosome}, Location: {chromosome_location}")
    else:
        gene_id_list.append("Not found")
        chromosome, chromosome_location, gene_location, gene_summary = "Not available", "Not available", "Not available", "Not available"
        print(f"Gene: {gene_name} not found")
    
    chromosome_list.append(chromosome)
    chromosome_location_list.append(chromosome_location)
    gene_location_list.append(gene_location)
    gene_summary_list.append(gene_summary)
    
    # 添加延迟以避免过快请求
    time.sleep(1)

# 将信息添加到DataFrame
df['GeneID'] = gene_id_list
df['Chromosome'] = chromosome_list
df['Chromosome Location'] = chromosome_location_list
df['Gene Location'] = gene_location_list
df['Gene Summary'] = gene_summary_list

# 保存更新后的DataFrame到新的Excel文件
output_file_path = r'e:/desk/基因注释_更新_完整信息.xlsx'
df.to_excel(output_file_path, index=False)

print(f"Process completed. Results saved to {output_file_path}")