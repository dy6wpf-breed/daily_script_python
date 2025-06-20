from Bio import Entrez
import pandas as pd

# Set the email for NCBI Entrez
Entrez.email = "bjlmwpf@163.com"

# Function to fetch detailed annotation from NCBI for a given GeneID
def fetch_gene_annotation(gene_id):
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        # Extracting detailed annotation from the fetched records
        if records:
            gene_info = records[0]
            # Extract the gene description
            gene_description = gene_info.get('Entrezgene_gene', {}).get('description', "No description available")
            return gene_description
        else:
            return "No gene information available"
    except Exception as e:
        return str(e)

# Load the Excel file
file_path = r'E:\desk\基因注释1.xlsx'  # Update this path
df = pd.read_excel(file_path)

# Fetch annotations for each GeneID
gene_ids = df['NCBI GeneID'].tolist()
detailed_annotations = [fetch_gene_annotation(str(gene_id)) for gene_id in gene_ids]

# Add the fetched annotations to the dataframe
df['Gene Description'] = detailed_annotations

# Save the updated dataframe to a new Excel file
output_file_path = r'E:\desk\基因注释1_更新_基因描述.xlsx'  # Update this path
df.to_excel(output_file_path, index=False)

print(f"Updated file saved to: {output_file_path}")