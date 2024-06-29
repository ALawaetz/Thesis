import os
import pandas as pd
import sys

the_ratio = float(sys.argv[1])
the_cutoff = float(sys.argv[2])

output_folder = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios'


file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/1d_pyReadCounters_RPKM_normalised_to_RNAseq/Jag_merged.novo_count_output_cDNAs.gtf_hittable_cDNAs.txt'
jag = pd.read_csv(file, sep = '\t', header = 0, usecols = ['Gene', 'seq_normalised_count'])

file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/1d_pyReadCounters_RPKM_normalised_to_RNAseq/KhpA_merged.novo_count_output_cDNAs.gtf_hittable_cDNAs.txt'
khpA = pd.read_csv(file, sep = '\t', header = 0, usecols = ['Gene', 'seq_normalised_count'])

file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/1d_pyReadCounters_RPKM_normalised_to_RNAseq/Kre_merged.novo_count_output_cDNAs.gtf_hittable_cDNAs.txt'
kre = pd.read_csv(file, sep = '\t', header = 0, usecols = ['Gene', 'seq_normalised_count'])

# Merge dataframes
merged_df = jag.merge(khpA, on='Gene', how='outer', suffixes=('_jag', '_khpA')).merge(kre, on='Gene', how='outer')
merged_df.columns = ['Gene', 'RPKMPW_jag', 'RPKMPW_khpA', 'RPKMPW_kre']

### Fill NaN values with 0
merged_df.fillna(0, inplace=True)

### for each HTF gene and each target, calculate HTF(gene)/min(HTF(gene))
jag_ratios = []
khpA_ratios = []
kre_ratios = []
for index, row in merged_df.iterrows():
    try:
        jag_ratios.append(row['RPKMPW_jag'] / min([row['RPKMPW_jag'], row['RPKMPW_khpA'], row['RPKMPW_kre']]))
    except ZeroDivisionError:
        jag_ratios.append(the_ratio)
    try:
        khpA_ratios.append(row['RPKMPW_khpA'] / min([row['RPKMPW_jag'], row['RPKMPW_khpA'], row['RPKMPW_kre']]))
    except ZeroDivisionError:
        khpA_ratios.append(the_ratio)
    try:
        kre_ratios.append(row['RPKMPW_kre'] / min([row['RPKMPW_jag'], row['RPKMPW_khpA'], row['RPKMPW_kre']]))
    except ZeroDivisionError:
        kre_ratios.append(the_ratio)

merged_df['jag_ratio'] = jag_ratios
merged_df['khpA_ratio'] = khpA_ratios
merged_df['kre_ratio'] = kre_ratios

### Sort dataframes by ratio cutoffs and RPKMPW cutoffs
merged_df[(merged_df['jag_ratio'] >= the_ratio) & (merged_df['RPKMPW_jag'] > the_cutoff)].sort_values(['RPKMPW_jag'], ascending = False).to_csv(os.path.join(output_folder, 'Jag_ratios.txt'), sep = '\t', header = True, index = False)
merged_df[(merged_df['khpA_ratio'] >= the_ratio) & (merged_df['RPKMPW_khpA'] > the_cutoff)].sort_values(['RPKMPW_khpA'], ascending = False).to_csv(os.path.join(output_folder, 'KhpA_ratios.txt'), sep = '\t', header = True, index = False)
merged_df[(merged_df['kre_ratio'] >= the_ratio) & (merged_df['RPKMPW_kre'] > the_cutoff)].sort_values(['RPKMPW_kre'], ascending = False).to_csv(os.path.join(output_folder, 'Kre_ratios.txt'), sep = '\t', header = True, index = False)




### One idea for ordering was calculating z_scores but this idea was trown as it was more simple
### to just calculate ratios between RPKMPW values as done above

# # Calculate row-wise mean and standard deviation
# merged_df['row_mean'] = merged_df[['RPKMPW_jag', 'RPKMPW_khpA', 'RPKMPW_kre']].mean(axis=1)
# merged_df['row_std'] = merged_df[['RPKMPW_jag', 'RPKMPW_khpA', 'RPKMPW_kre']].std(axis=1)
# # Calculate z-scores for each protein
# merged_df['Z_Score_jag'] = (merged_df['RPKMPW_jag'] - merged_df['row_mean']) / merged_df['row_std']
# merged_df['Z_Score_khpA'] = (merged_df['RPKMPW_khpA'] - merged_df['row_mean']) / merged_df['row_std']
# merged_df['Z_Score_kre'] = (merged_df['RPKMPW_kre'] - merged_df['row_mean']) / merged_df['row_std']
# # Drop the temporary columns used for mean and std deviation
# merged_df.drop(columns=['row_mean', 'row_std'], inplace=True)


