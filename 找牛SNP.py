from Bio import Entrez
import pandas as pd
import os

# NCBI Entrez 设置
Entrez.email = "dy6wpf@gmail.com"  # 替换为你的电子邮件

def fetch_gene_info(gene_id):
    try:
        # 获取基因信息
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        # 调试信息
        print(f"Fetched data for gene ID {gene_id}: {records}")

        # 解析基因信息
        if records and 'Entrezgene' in records:
            gene_info = records['Entrezgene'][0]
            gene_name = gene_info.get('Entrezgene_gene', {}).get('Gene-ref', {}).get('Gene-ref_locus', [''])[0]
            return gene_name
        else:
            print(f"No data found for gene ID {gene_id}")
            return None
    except Exception as e:
        print(f"Error fetching data for {gene_id}: {e}")
        return None

def main():
    genes = [
        "ENSBTAG00000054211", "DNAH11", "CDCA7L", "IL6", "TOMM7", "CCDC126",
        "DBF4", "SLC25A40", "ELAPOR2", "GRM3", "RELN", "PIK3CG", "NRCAM",
        "ST7", "GPR85", "TMEM168", "DOCK4", "IMMP2L", "TPK1", "KRR1",
        "OSBPL8", "SORCS2", "FGFR3", "NMUR2", "TXNDC8", "ENSBTAG00000052141",
        "METTL24", "TMEM182", "SH3RF3", "CSMD3", "PIP4P2", "CNGB3",
        "ENSBTAG00000053138", "CD82", "TP53I11", "PHF21A", "AMBRA1",
        "ENSBTAG00000048943", "TMIGD1", "CPD", "GOSR1", "ABCC3", "ATP6V0A1",
        "ADAMTS6", "MITF", "F13A1", "PRDM5", "JAG1", "SEC61G", "GCLM",
        "ENSBTAG00000026163", "ENSBTAG00000023541", "ENSBTAG00000031913",
        "ENSBTAG00000038619", "ENSBTAG00000054689", "ENSBTAG00000054580",
        "ENSBTAG00000047461", "ENSBTAG00000031458", "ENSBTAG00000053188",
        "ENSBTAG00000048304", "RIOK2", "MEIS1", "ANTXR1", "PRKG1",
        "UBE2E3", "CLCA2", "LARGE1", "DNM1L", "CCDC91", "SOX5", "ETV6",
        "TAS2R42"
    ]

    results = []
    for gene_id in genes:
        gene_name = fetch_gene_info(gene_id)
        if gene_name:
            results.append({"Gene": gene_id, "Gene_Name": gene_name})

    # 如果没有结果，输出提示
    if not results:
        print("No gene information was retrieved.")

    # 保存结果到指定目录
    output_directory = 'E:\\desk'
    output_file = os.path.join(output_directory, 'gene_info_results.csv')
    
    # 创建目录（如果不存在）
    os.makedirs(output_directory, exist_ok=True)
    
    # 将结果保存到 CSV 文件
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
