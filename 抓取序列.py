from Bio import Entrez
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import argparse

# 设置数据库列表，用户可以选择搜索的数据库
databases = [
    'pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest',
    'structure', 'genome', 'books', 'cancerchromosomes', 'cdd', 'gap',
    'domains', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
    'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pmc',
    'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound',
    'pcsubstance', 'snp', 'taxonomy', 'toolkit', 'unigene', 'unists'
]

# 解析命令行参数
parser = argparse.ArgumentParser(description='This script is used to fetch sequences from NCBI.')
parser.add_argument('-t', '--term', help='Input search term, e.g., "cow"', required=True)
parser.add_argument('-d', '--database', help='Database to search, default is nucleotide', default='nucleotide', required=False)
parser.add_argument('-r', '--rettype', help='Return type, fasta or gb, default is gb', default='gb', required=False)
parser.add_argument('-o', '--out_dir', help='Output directory path, default is current working directory', default=os.getcwd(), required=False)
parser.add_argument('-n', '--name', help='Specify the output file name, default is "seq"', default='seq', required=False)

args = parser.parse_args()

# 检查输出目录是否存在，如果不存在则创建
if os.path.exists(args.out_dir):
    output_dir = os.path.abspath(args.out_dir)
else:
    os.makedirs(args.out_dir)
    output_dir = os.path.abspath(args.out_dir)

# 设置输出文件的完整路径
output_file_path = os.path.join(output_dir, f"{args.name}.{args.rettype}")

# 打开输出文件
with open(output_file_path, "w") as output_handle:
    # 设置NCBI的电子邮件地址
    Entrez.email = "bjlmwpf@163.com"  # 替换为你的电子邮件地址

    # 执行搜索
    handle = Entrez.esearch(db=args.database, term=args.term, idtype="acc")
    record = Entrez.read(handle)
    handle.close()

    # 遍历搜索结果中的ID列表
    for id in record['IdList']:
        print(f"Fetching sequence for ID: {id}")
        # 获取序列数据
        handle = Entrez.efetch(db=args.database, id=id, rettype=args.rettype, retmode="text")
        record = SeqIO.read(handle, args.rettype)
        handle.close()
        # 写入到输出文件
        SeqIO.write(record, output_handle, args.rettype)

print(f"Sequences have been fetched and saved to {output_file_path}")