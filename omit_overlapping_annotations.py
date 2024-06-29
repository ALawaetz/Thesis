import pandas as pd
import sys

input_gtf = sys.argv[1]
output_gtf = f'{input_gtf[:-4]}_no_overlaps.gtf'


def process_gtf(input_gtf, feature):
    # Read the GTF file into a DataFrame
    df = input_gtf

    # Filter rows with feature 'sRNA'
    srna_df = df[df[1] == feature]

    # Identify overlapping 'sRNA' annotations
    overlaps = []
    for idx1, row1 in srna_df.iterrows():
        for idx2, row2 in srna_df.iterrows():
            if idx1 != idx2 and overlap(row1, row2):
                overlaps.append((idx1, idx2))

    # Keep only the longest annotation for each overlap
    to_drop = set()
    for idx1, idx2 in overlaps:
        if idx1 not in to_drop and idx2 not in to_drop:
            if row_length(srna_df.loc[idx1]) > row_length(srna_df.loc[idx2]):
                to_drop.add(idx2)
            else:
                to_drop.add(idx1)

    # Drop the rows to be removed
    df_filtered = df.drop(to_drop)

    # Save the modified GTF file
    return df_filtered
    

def overlap(row1, row2):
    # Check if two GTF rows overlap
    return (row1[3] <= row2[4] and row1[4] >= row2[3] and row1[6] == row2[6]) or (row2[3] <= row1[4] and row2[4] >= row1[3] and row1[6] == row2[6])

def row_length(row):
    # Calculate the length of a GTF row
    return row[4] - row[3]

# Example usage:
input_gtf = pd.read_csv(input_gtf, sep='\t', header=None)
features = list(set(input_gtf[1].tolist()))
for feature in features:
    if feature != 'CDS':
        print(feature)
        input_gtf = process_gtf(input_gtf, feature)
    else:
        pass

input_gtf.to_csv(output_gtf, sep='\t', header=False, index=False)
