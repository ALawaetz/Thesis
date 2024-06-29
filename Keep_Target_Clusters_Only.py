import os
import pandas as pd
import gff3
import sys

cluster_file = sys.argv[1]
cluster_df = pd.read_csv(cluster_file, sep = '\t', names = gff3.header, comment = '#')

ratio_file = sys.argv[2]
ratio_df = pd.read_csv(ratio_file, sep = '\t', header = 0)

output_folder = sys.argv[3]


for index, row in cluster_df.iterrows():
    attribute_list = row['attributes'].split(';')
    for item in attribute_list:
        if 'gene_name' in item:
            gene_name_list = item.split('gene_name ')[1].split(',')
            break
        else:
            pass
    gene_name_list = [i.strip('"') for i in gene_name_list]
    x = 0
    for gene in gene_name_list:
        if gene in ratio_df['Gene'].tolist():
            x += 1
    if x != len(gene_name_list):
        cluster_df = cluster_df[cluster_df.index != index]

cluster_df.to_csv(os.path.join(output_folder, os.path.basename(cluster_file)), sep = '\t', header = None, index = False)