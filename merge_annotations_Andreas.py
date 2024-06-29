import os
import sys
import pandas as pd
import gff3
import re

### This script takes a gtf file as input and merges overlapping annotations

file = sys.argv[1]

df = pd.read_csv(file, sep = '\t', comment = '#', names = gff3.header)
df = df.sort_values(['start', 'end'], ascending = True)
df_fwd = df[df['strand'] == '+'].reset_index(drop = True)
df_rev = df[df['strand'] == '-'].reset_index(drop = True)

merged_df = pd.DataFrame()

def merge_func(df, the_strand):
    start = []
    end = []
    attributes = []
    n = 0
    while n <= max(df.index):
        start.append(df['start'].tolist()[n])
        att_temp = df['attributes'].tolist()[n]
        try:
            while df['start'].tolist()[n + 1] <= df['end'].tolist()[n]:
                n += 1
                att_temp += df['attributes'].tolist()[n]
            end.append(df['end'].tolist()[n])
        except IndexError:
            end.append(df['end'].tolist()[n])
        n += 1
        att_temp2 = ''
        for att in att_temp.split(';'):
            if 'gene_name' in att:
                att_temp2 += att.split('gene_name ')[1].strip('"')
                att_temp2 += ','
        att_temp2 = att_temp2.rstrip(',')
        alist = list(set(att_temp2.split(',')))
        alist = ','.join(alist)
        attributes.append(f'gene_name {alist}')

    merged_df = pd.DataFrame()
    merged_df['id'] = ['NC_000964.3'] * len(start)
    merged_df['feature'] = ['cluster'] * len(start)
    merged_df['source'] = ['interval'] * len(start)
    merged_df['start'] = start
    merged_df['end'] = end
    merged_df['max_height'] = ['.'] * len(start)
    merged_df['strand'] = [the_strand] * len(start)
    merged_df['frame'] = ['.'] * len(start)
    merged_df['attributes'] = attributes
    return merged_df

merge_fwd = merge_func(df_fwd, '+')
merge_rev = merge_func(df_rev, '-')
final_merged_df = pd.concat([merge_fwd, merge_rev])

final_merged_df.to_csv(file + 'merged', sep = '\t', header = None, index = False)


### This code adds "" around gene names
def process_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            # Use regex to find the gene_name field and surround it with quotes
            modified_line = re.sub(r'(gene_name\s)(.+)', r'\1"\2"', line)
            outfile.write(modified_line)

# Example usage
input_filename = file + 'merged'
output_filename = file + "final_merge"
process_file(input_filename, output_filename)

