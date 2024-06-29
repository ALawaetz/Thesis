import os
import pandas as pd
import sys
# import pool_annotations_within_pyRAP_CRAC_modified


file = '/Users/andreas/Bacillus/Bioinformatics/pyRAP-main-2/pyRAP_output_B_subtilis/Pipeline/Annotation_files/All_conditions.gff3_categorized_operons.gff3_supplied_with_SubtiWiki_Genbank_w_predicted_UTRs.gff3.gff3'
df = pd.read_csv(file, sep = '\t', header = None)


### Change 5'UTR to 5UTR and 3'UTR to 3UTR
features = []

for feature, source in zip(df[2].tolist(), df[1].tolist()):
    if feature == "5'UTR":
        features.append('5UTR')
    elif feature == "3'UTR":
        features.append('3UTR')
    else:
        features.append(feature)


df[2] = features


### Change RendSeq sRNA features to sRNA_type e.g. sRNA_independent
sRNA_type = []
for feature, attribute in zip(df[2].tolist(), df[8].tolist()):
    if feature != 'sRNA':
        sRNA_type.append('.')
    else:
        alist = attribute.split(';')
        z = 0
        for a in alist:
            if a.find(' ') == 0:
                a = a[1:]
            else:
                pass
            if a.find('sRNA_type') == 0:
                sRNA_type.append(a[10:])
                z += 1
                break
            else:
                pass
        if z == 0:
            print(alist)
            


features = []
if len(sRNA_type) == len(df):
    for feature, type in zip(df[2].tolist(), sRNA_type):
        if feature != 'sRNA':
            features.append(feature)
        else:
            if type == "5'UTR":
                features.append('{}_5UTR'.format(feature))
            elif type == "3'UTR":
                features.append('{}_3UTR'.format(feature))
            elif type == "independent_w_isoform":
                features.append('{}_independent'.format(feature))
            elif type == "independent":
                features.append('{}_independent'.format(feature))
            elif type == "inter/intragenic":
                features.append('{}_inter_intragenic'.format(feature))
            elif type == "intergenic":
                features.append('{}_inter_intragenic'.format(feature))
            elif type == "intragenic":
                features.append('{}_inter_intragenic'.format(feature))
            elif type == 'non_sRNA/non_CDS_overlap':
                features.append('{}_other'.format(feature))
            else:
                print('review code')
                print(type)
                sys.exit()
else:
    print('#######review code')

df[2] = features



### Move features to Source column
df[1] = df[2].tolist()

### Change feature column so only uses CDS and exon
df[2] = df[2].apply(lambda x: 'CDS' if x == 'CDS' else 'exon')


### Redo pool module from pyRAP in a new simplified version
df = df.rename(columns={0: 'seqID', 1: 'source', 2: 'feature', 3: 'start', 4: 'end', 5: 'score', 6: 'strand', 7: 'phase', 8: 'attributes'})

def remove_duplicates(df):
    # Sort the DataFrame by Column 1, Column 2, and Column 3
    df_sorted = df.sort_values(by=['source', 'start', 'end']).reset_index(drop=True)

    # Initialize an empty list to store the indices of rows to be removed
    rows_to_remove = []

    # Iterate over the sorted DataFrame
    for i in range(1, len(df_sorted)):
        prev_row = df_sorted.iloc[i - 1]
        curr_row = df_sorted.iloc[i]

        # Check if Column 1 is the same and Column 2 and Column 3 are within +/- 5
        if (
            prev_row['source'] == curr_row['source'] and
            abs(prev_row['start'] - curr_row['start']) <= 5 and
            abs(prev_row['end'] - curr_row['end']) <= 5
        ):
            rows_to_remove.append(i)

    # Drop the rows marked for removal
    df_modified = df_sorted.drop(rows_to_remove)
    df_discarded = df_sorted[df_sorted.index.isin(rows_to_remove)]

    return [df_modified, df_discarded]

df_modified = remove_duplicates(df)[0]

df_modified = df_modified.sort_values(by = ['start', 'end'])



### Change attributes column so only contains gene_name and gene_id
### gene_id shall be LocusTag if possible or else genomic coordinates

locus = []
names = []
for item in df_modified['attributes'].tolist():
    alist = item.split(';')
    x = 0
    for a in alist:
        if a.find('LocusTag') == 0: #beware that some files are spelled LocusTag with capital L and other files with small l
            locus.append(a[9:])
            x += 1
            break
        else:
            pass
    if x == 0:
        locus.append('no_locustag')
    x = 0
    for a in alist:
        if a.find('Name=') == 0:
            names.append(a[5:])
            x += 1
            break
        else:
            pass
    if x == 0:
        for a in alist:
            if a.find('name=') == 0:
                names.append(a[5:])
                x += 1
                break
            else:
                pass
    if x == 0:
        for a in alist:
            if a.find('ID=') == 0:
                names.append(a[3:])
                x += 1
                break
            else:
                pass
    if x == 0:
        names.append('no_name')

### Neither pyCRAC or DESeq2 does well with names containing spaces, so change those
names = [i.replace(' ', '@') for i in names]

### Redo names, so no names are the same. Hence, e.g. dnaA.1, dnaA.2 and so forth.
newnames = []
numbers = []
for item in names:
    if item in newnames:
        number = newnames.count(item)
        newnames.append(item)
        numbers.append(number + 1)
    else:
        newnames.append(item)
        numbers.append(1)

version = []
for new, num in zip(newnames, numbers):
    version.append('{}.{}'.format(new, num))


attributes = []
for name, locus in zip(version, locus):
    attributes.append('gene_name={}; locusTag={}'.format(name, locus))

df_modified['attributes'] = attributes



### Finally, it just needs to be converted to gtf to run the CRAC_pipeline_PE.py

df_modified['attributes'] = [i.replace('=', ' ') for i in df_modified['attributes'].tolist()]

df_modified.to_csv(file + '_CRAC_modified.gtf', sep = '\t', header = None, index = False)





