import os
import glob
import pandas as pd
import sys




### reads file
reads_file = sys.argv[1]

### cluster file
cluster_file = sys.argv[2]
cluster_df = pd.read_csv(cluster_file, sep = '\t', header = None, comment = '#')

gtf_file = sys.argv[3]

root_folder = sys.argv[4]

output_folder = sys.argv[5]

feature = sys.argv[6]

ratios_file = sys.argv[7]

top_genes = int(float(sys.argv[8]))

# ratios = pd.read_csv(ratios_file, sep = '\t', usecols = ['Gene'], header = 0)
# include_genes = ratios['Gene'].tolist()[:top_genes]

include_genes = []
for index, row in cluster_df.iterrows():
    attribute_list = row[8].split(';')
    for item in attribute_list:
        if 'gene_name' in item:
            gene_name_list = item.split('gene_name ')[1].split(',')
            break
        else:
            pass
    gene_name_list = [i.strip('"') for i in gene_name_list]
    include_genes += gene_name_list



### make GTF file only containing genes in top100 list
### extract relevant features
pyRAP = pd.read_csv(gtf_file, sep = '\t', header = None)
names = []
for item in pyRAP[8].tolist():
    alist = item.split(' ')
    names.append(alist[1])
pyRAP['name'] = names
pyRAP = pyRAP[pyRAP['name'].isin(include_genes)]
pyRAP.to_csv(os.path.join(root_folder, f'{os.path.basename(reads_file)[:-4]}_temp.gtf'), sep = '\t', header = None, index = False)
gtf_file = root_folder + f'/{os.path.basename(reads_file)[:-4]}_temp.gtf'


if feature == 'All':
    os.system('python pyBinCollector.py -f {} --gtf {} -n 50 --normalize --deletions -o {}/{}'.format(reads_file, gtf_file, output_folder, os.path.basename(reads_file)))
elif feature == '5UTR' or feature == '3UTR':
    os.system('python pyBinCollector.py -f {} --gtf {} -a {} -n 20 --normalize --deletions -o {}/{}'.format(reads_file, gtf_file, feature, output_folder, os.path.basename(reads_file)))
else:
    os.system('python pyBinCollector.py -f {} --gtf {} -a {} -n 100 --normalize --deletions -o {}/{}'.format(reads_file, gtf_file, feature, output_folder, os.path.basename(reads_file)))

os.remove(gtf_file)


