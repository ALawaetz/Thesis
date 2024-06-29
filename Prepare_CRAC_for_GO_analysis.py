import os
import pandas as pd
import glob


### prepare CRAC to GO analysis
file = '/Users/andreas/Bacillus/Bioinformatics/Scripts_RNAseq/subtiwiki.gene.export.2024-05-25.tsv'
subti = pd.read_csv(file, sep = '\t', header = 0)


file = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files/B_subtilis_GenBank_w_name_and_locus.gff3'
genbank = pd.read_csv(file, sep = '\t', header = 0, comment = '#')

### some genes are on the form geneX_geneY
### split these genes into two rows; geneX and geneY

input_folder = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios'
output_folder = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/ratios/prepared'


for file in glob.glob(input_folder + '/*.txt'):
    df = pd.read_csv(file, sep = '\t', header = 0)
    gene_base = []
    for gene in df['Gene'].tolist():
        if gene.find('NP_') == 0 and gene.count('_') == 1:
            gene_base.append(gene.split('.')[0] + '.' + gene.split('.')[1])
        elif gene.find('YP_') == 0 and gene.count('_') == 1:
            gene_base.append(gene.split('.')[0] + '.' + gene.split('.')[1])
        else:
            gene_base.append(gene.split('.')[0])
    df['Gene_base'] = gene_base
    df = df.drop_duplicates(['Gene_base'])
    df.to_csv(file, sep = '\t', header = True, index = False)

for file in glob.glob(input_folder + '/*.txt'):
    all = pd.read_csv(file, sep = '\t', header = 0)

    ### replace 'BSU_' with 'BSU^' 
    ### This way we can seperate transcripts by _ followingly
    temp_gene_names = []
    for name in all['Gene_base'].tolist():
        temp_gene_names.append(name.replace('BSU_', 'BSU^').replace('NP_', 'NP^').replace('YP_', 'YP^').replace('Novel_transcript', 'Novel^transcript'))
    all['Gene_base'] = temp_gene_names

    _df = pd.DataFrame()

    # Loop through each row in the original DataFrame
    for index, row in all.iterrows():
        # Check if 'Gene_base' contains an underscore
        if '_' in row['Gene_base']:
            # Split the 'Gene_base' value into a list of genes
            gene_list = row['Gene_base'].split('_')
            # Create new rows for each gene in gene_list
            for gene in gene_list:
                new_row = pd.DataFrame({
                    'Gene_base': [gene]
                })
                _df = pd.concat([_df, new_row], ignore_index=True)
        ### replace '@' with ' ' in gene names
        elif '@' in row['Gene_base']:
            new_row = pd.DataFrame({
                    'Gene_base': [row['Gene_base'].replace('@', ' ')]
                })
            _df = pd.concat([_df, new_row], ignore_index=True)

        else:
            # If nothing else, keep the row as is
            all_temp = pd.DataFrame()
            all_temp['Gene_base'] = all['Gene_base'][all.index == index].tolist()
            _df = pd.concat([_df, all_temp], ignore_index=True)
    
    
    ### replace 'BSU_' with 'BSU^' 
    ### This way we can seperate transcripts by _ followingly
    revert_gene_names = []
    for name in _df['Gene_base'].tolist():
        revert_gene_names.append(name.replace('BSU^', 'BSU_').replace('NP^', 'NP_').replace('YP^', 'YP_').replace('Novel^transcript', 'Novel_transcript'))
    _df['Gene_base'] = revert_gene_names
    _df = _df.drop_duplicates(['Gene_base'])

    ### make list of DEG
    sig_genes = _df['Gene_base'].tolist()

    ### convert gene names to BSU_numbers
    bsu_numbers = []
    # for gene in sig_genes:
    #     try:
    #         bsu_numbers.append(subti['locus'][subti['title'] == gene].iloc[0])
    #     except IndexError:
    #         x = 0
    #         for index, row in subti.iterrows():
    #             try:
    #                 if gene in row['synonyms']:
    #                     bsu_numbers.append(subti['locus'][subti.index == index].iloc[0])
    #                     x += 1
    #                     break
    #             ## for na values
    #             except TypeError:
    #                 pass
    #         if x == 0:
    #             print(f'Error: Couldnt convert {gene}')
    for gene in sig_genes:
        try:
            bsu_numbers.append(subti['locus'][subti['title'] == gene].iloc[0])
            continue
        except IndexError:
            x = 0
            for index, row in subti.iterrows():
                try:
                    if gene in row['synonyms']:
                        bsu_numbers.append(subti['locus'][subti.index == index].iloc[0])
                        x += 1
                        break
                ## for na values
                except TypeError:
                    pass
            if x == 0:
                try:
                    bsu_numbers.append(genbank['locusTag'][genbank['name'] == gene].iloc[0])
                    continue
                except IndexError:
                    print(f'Error: Couldnt convert {gene}')
                    bsu_numbers.append('no_locusTag')
    
    bsunumbers = [i.replace('_', '') for i in bsu_numbers]

    ### convert bsu_numbers to genbank names
    genbank_names = []
    for bsu in bsu_numbers:
        for index, att in enumerate(genbank['attributes'].tolist()):
            x = 0
            if bsu in att:
                att_list = att.split(';')
                for a in att_list:
                    if a.find('gene=') == 0:
                        genbank_names.append(a[5:])
                        x += 1
                        break
            if x == 1:
                break
    
    ### Write list to .txt file
    # Specify the file name
    file_name = f"{output_folder}/{os.path.basename(file)}_sig_genes.txt"

    # Open the file in write mode
    with open(file_name, "w") as filex:
        # Iterate through each item in the list
        for item in genbank_names:
            # Write the item followed by a newline character
            filex.write(f"{item}\n")

    file_name = f"{output_folder}/{os.path.basename(file)}_sig_genes_BSU_Numbers.txt"

    # Open the file in write mode
    with open(file_name, "w") as filex:
        # Iterate through each item in the list
        for item in bsu_numbers:
            # Write the item followed by a newline character
            filex.write(f"{item}\n")

    file_name = f"{output_folder}/{os.path.basename(file)}_sig_genes_BSUNumbers_no_underscore.txt"

    # Open the file in write mode
    with open(file_name, "w") as filex:
        # Iterate through each item in the list
        for item in bsunumbers:
            # Write the item followed by a newline character
            filex.write(f"{item}\n")
