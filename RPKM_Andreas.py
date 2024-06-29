import os
import pandas as pd
import sys
import gff3

### Process pyReadCounters file
### Calculate RPKM

pyReadCounter_file = sys.argv[1]
gtf_file = sys.argv[2]
output_folder = sys.argv[3]
coverageBed_file = sys.argv[4]
cluster_height = int(float(sys.argv[5]))

df = pd.read_csv(pyReadCounter_file, sep = '\t', names = ['Gene', 'sense', 'antisense', 'cDNA_count'], comment = '#')
df = df[df['sense'] >= cluster_height]
gtf = pd.read_csv(gtf_file, sep = '\t', names = gff3.header)
wt = pd.read_csv(coverageBed_file, sep = '\t', names = ['Gene', 'Count'])

### Calculate RPKM
genes = []
for attribute in gtf['attributes'].tolist():
    genes.append(attribute[10:])
gtf['Gene'] = genes
RPKM = []
depth = sum(df['cDNA_count'])
for index, row in df.iterrows():
    gene_length = gtf['end'][gtf['Gene'] == row['Gene']].iloc[0] - gtf['start'][gtf['Gene'] == row['Gene']].iloc[0]
    RPKM.append(int(row['cDNA_count'] * 1000 * 1000000 / (gene_length * depth)))
df['RPKM'] = RPKM
df = df.drop(['sense', 'antisense', 'cDNA_count'], axis = 1)


### normalised RPKM to WT RNA-seq
normalised_count = []
for index, row in df.iterrows():
    count_value = wt['Count'][wt['Gene'] == row['Gene']].iloc[0]
    if count_value >= 100:
        norm = row['RPKM'] / count_value
    else:
        norm = 0
    normalised_count.append(norm)

df['seq_normalised_count'] = normalised_count

df.to_csv(os.path.join(output_folder, os.path.basename(pyReadCounter_file)), sep = '\t', header = True, index = False)
