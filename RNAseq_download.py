import requests
import csv
import subprocess
import os

def get_gene_id(protein_id):
    """
    使用Ensembl的xrefs API将蛋白质ID映射为基因ID
    """
    url = f"https://rest.ensembl.org/xrefs/symbol/selaginella_moellendorffii/{protein_id}?"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        mappings = response.json()
        if mappings:
            # 找到相应的基因ID
            return mappings[0]['id']
        else:
            print(f"No mappings found for {protein_id}")
            return None
    else:
        print(f"Failed to retrieve mappings: HTTP {response.status_code}")
        return None

def get_gene_location(gene_id):
    """
    从Ensembl获取基因的位置信息
    """
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        gene_info = response.json()
        chromosome = gene_info['seq_region_name']
        start = gene_info['start']
        end = gene_info['end']
        strand = 'positive' if gene_info['strand'] == 1 else 'negative'
        exons = gene_info.get('Exon', [])
        exon_positions = [(exon['start'], exon['end']) for exon in exons]
        return chromosome, start, end, strand, exon_positions
    else:
        print(f"Failed to retrieve gene location: HTTP {response.status_code}")
        return None, None, None, None, None

def fetch_sequence(chromosome, start, end):
    """
    根据位置从Ensembl提取基因组序列
    """
    url = f"https://rest.ensembl.org/sequence/region/selaginella_moellendorffii/{chromosome}:{start}..{end}:1?"
    headers = {"Content-Type": "text/plain"}
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to fetch sequence: HTTP {response.status_code}")
        return None

def grep_gene_info(gene_id, gff_file):
    """
    使用grep从GFF文件中提取特定基因的信息
    """
    command = f"grep {gene_id} {gff_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout

def parse_gff_output(gff_output):
    """
    解析grep命令提取的GFF信息，并将exon、CDS、UTR和mRNA信息组合起来
    """
    gff_lines = gff_output.strip().split("\n")
    cds = []
    exons = []
    five_prime_utrs = []
    three_prime_utrs = []
    mrna = []

    for line in gff_lines:
        columns = line.split("\t")
        feature_type = columns[2]
        start = columns[3]
        end = columns[4]

        if feature_type == "exon":
            exons.append([start, end])
        elif feature_type == "CDS":
            cds.append([start, end])
        elif feature_type == "five_prime_UTR":
            five_prime_utrs.append([start, end])
        elif feature_type == "three_prime_UTR":
            three_prime_utrs.append([start, end])
        elif feature_type == "mRNA":
            mrna.append([start, end])

    # 去重处理
    cds = list(set(tuple(interval) for interval in cds))
    exons = list(set(tuple(interval) for interval in exons))
    five_prime_utrs = list(set(tuple(interval) for interval in five_prime_utrs))
    three_prime_utrs = list(set(tuple(interval) for interval in three_prime_utrs))
    mrna = list(set(tuple(interval) for interval in mrna))

    return cds, exons, five_prime_utrs, three_prime_utrs, mrna

def process_gene_list(input_file, output_file, gff_file):
    """
    处理基因列表，将从Ensembl提取的序列以及从GFF文件中提取的相关信息写入CSV文件
    """
    with open(input_file, 'r') as file:
        gene_ids = file.read().splitlines()

    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Gene ID', 'Sequence', 'Chromosome', 'Strand', 'Start', 'End', 'Exon Count', 'Exon Positions',
                      'CDS', 'Exons', 'mRNA', '5\' UTRs', '3\' UTRs']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for gene_id in gene_ids:
            print(f"Processing {gene_id}...")

            # 先尝试将蛋白质ID映射为基因ID（如果需要）
            mapped_gene_id = get_gene_id(gene_id)
            if mapped_gene_id:
                gene_id = mapped_gene_id

            # 从Ensembl提取基因位置信息
            chromosome, start, end, strand, ensembl_exons = get_gene_location(gene_id)
            if chromosome and start and end:
                extended_start = start - 5000
                extended_end = end + 5000
                sequence = fetch_sequence(chromosome, extended_start, extended_end)

                exon_count = len(ensembl_exons)
                exon_positions = "; ".join([f"{start}-{end}" for start, end in ensembl_exons])

                # 从GFF文件提取CDS、exon、UTR、mRNA信息
                gff_output = grep_gene_info(gene_id, gff_file)
                cds, exons, five_prime_utrs, three_prime_utrs, mrna = parse_gff_output(gff_output)

                row_data = {
                    'Gene ID': gene_id,
                    'Sequence': sequence,
                    'Chromosome': chromosome,
                    'Strand': strand,
                    'Start': start,
                    'End': end,
                    'Exon Count': exon_count,
                    'Exon Positions': exon_positions,
                    'CDS': ' | '.join([f"{start}:{end}" for start, end in cds]),
                    'Exons': ' | '.join([f"{start}:{end}" for start, end in exons]),
                    'mRNA': ' | '.join([f"{start}:{end}" for start, end in mrna]),
                    '5\' UTRs': ' | '.join([f"{start}:{end}" for start, end in five_prime_utrs]),
                    '3\' UTRs': ' | '.join([f"{start}:{end}" for start, end in three_prime_utrs])
                }

                writer.writerow(row_data)
            else:
                print(f"Skipping {gene_id} due to missing location information.")

# 设置输入文件和输出文件的路径
os.chdir('/Users/wangyuanyuan/Downloads/1028以后文件/database_poged/')
input_file = 'sm_gene.txt'
output_file = 'sm_gene_info.csv'
gff_file = 'Selaginella_moellendorffii.v1.0.59.gff3'

# 处理基因列表并写入CSV文件
process_gene_list(input_file, output_file, gff_file)

print("Finished processing all genes.")
