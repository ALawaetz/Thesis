import os
import glob
import shutil
import pandas as pd
import subprocess
import sys

### uncomment as you go along

raw_data_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/01.RawData/'
### all_in_one_folder should not be in the same directory as the raw_data_folder
all_in_one_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/All_in_one/'
reference_fasta = '/Users/andreas/Bacillus/Bioinformatics/FASTA/NC_000964v3.fasta'
analysis_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis'
sam_folder = f'{analysis_folder}/sam_files/'
bam_folder = f'{analysis_folder}/bam_files/'
reference_gff = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files/B_subtilis_CRAC_RNAseq_analysis.gff3'
CoverageBed_folder = f'{analysis_folder}/Coverage_bed'



# chrom_size = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/B_subtilis_chromsize.txt'
# base_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_spoVG_LB_OD3/Analysis/'
# sorted_bam_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_spoVG_LB_OD3/Analysis/Sorted_bam_files/'
# merged_bam_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_spoVG_LB_OD3/Analysis/merged_bam_folder'
# sorted_indexed_bam_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/sorted_indexed_bam_files'
# prepare_script = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/prepare_to_DESeq2.py'

all_folders = [analysis_folder, all_in_one_folder, sam_folder, bam_folder, CoverageBed_folder]

# Create folders if they don't already exist
for folder in all_folders:
    if not os.path.exists(folder):
        os.makedirs(folder)

### Move all files to the same folder
### Iterate over all folders in the raw_data_folder
# for root, dirs, files in os.walk(raw_data_folder):
#     for file in files:
#         if file.endswith('.fq.gz'):
#             # Construct full file path
#             file_path = os.path.join(root, file)
#             # Move the file to the target folder
#             shutil.move(file_path, all_in_one_folder)
#             print(f'Moved: {os.path.basename(file_path)} to {all_in_one_folder}')


### gunzip all gz files
# os.chdir(all_in_one_folder)
# os.system('gunzip *.gz')

### Build reference
# os.chdir(analysis_folder)
# os.system(f'bowtie2-build {reference_fasta} {os.path.basename(reference_fasta)}')

# ### reverse complement all reverse files
# for file in glob.glob(all_in_one_folder + '/*1.fq'):
#     os.system('seqtk seq -r {} > {}_flipped.fq'.format(file, file[:-3]))
#     os.system('rm {}'.format(file))


# ### align reads to reference genome using the build files
# file1 = []
# file2 = []
# for file in sorted(glob.glob(all_in_one_folder + '/*.fq')):
#     if file[-4] == '2':
#         file2.append(file)
#     else:
#         file1.append(file)
# for one, two in zip(file1, file2):
    # print(f'aligning {os.path.basename(file1)} and {os.path.basename(file2)} to genome')
#     os.system(f'bowtie2 -x {analysis_folder}/NC_000964v3.fasta -1 {two} -2 {one} -S {sam_folder}/{os.path.basename(two)[:-5]}.sam --ff --very-sensitive-local -p10')
    
### convert sam to bam
# for file in glob.glob(sam_folder + '/*.sam'):
#     print(f'converting {os.path.basename(file)} to bam format')
#     os.system('samtools view -bS {} > {}/{}.bam'.format(file, bam_folder, os.path.basename(file[:-4])))

# ## sort bam files
# for file in glob.glob(bam_folder + '/*.bam'):
#     print(f'sorting {os.path.basename(file)}')
#     os.system('samtools sort {} -T {} -o {}'.format(file, os.path.basename(file), file))


## merge bam files for samples too large to fit in one fq file
# bam_dict = {}
# for file in glob.glob(bam_folder + '/*.bam'):
#     if os.path.basename(file)[:7] in list(set(bam_dict.keys())):
#         bam_dict[os.path.basename(file)[:7]] += [file]
#     else:
#         bam_dict[os.path.basename(file)[:7]] = [file]

# for key, item in bam_dict.items():
#     if len(item) == 2:
#         print(f'merging {os.path.basename(item[0])} and {os.path.basename(item[1])}')
#         os.system('samtools merge {}_merged.bam {} {}'.format(item[0][:-4], item[0], item[1]))
#         os.system('rm {} {}'.format(item[0], item[1]))
#     elif len(item) > 2:
#         print('adjust code in samtools merge to accomodate more files')
#     else:
#         pass


### make coverageBed files to input DESeq2
for file in glob.glob(bam_folder + '/*.bam'):
    os.system('coverageBed -s -sorted -counts -a {} -b {} > {}/{}.txt'.format(reference_gff, file, CoverageBed_folder, os.path.basename(file)[:-4]))


# for file in glob.glob(CoverageBed_folder + '/*.txt'):
#     os.system('cut -f9-10 {} > {}_cut.txt'.format(file, file[:-4]))
#     os.system('rm {}'.format(file))


# # ### Only keep the name in the first column
# for file in glob.glob(CoverageBed_folder + '/*.txt'):
#     names = []
#     df = pd.read_csv(file, sep = '\t', header = None)
#     for item in df[0].tolist():
#         alist = item.split(';')
#         x = 0
#         for a in alist:
#             if a.find('Name=') == 0 or a.find('name=') == 0:
#                 names.append(a[5:])
#                 x += 1
#             else:
#                 pass
#         if x == 0:
#             print(item)
#             sys.exit()

#     df[0] = names

#     df.to_csv(file, sep = '\t', header = None, index = False)


# ### DESeq2 can't handle gene names containing space characters (' ')
# ### Unless it has quotation marks around
# for file_path in glob.glob(CoverageBed_folder + "/*.txt"):
#     subprocess.run(["bash", "-c", f'awk -F"\t" \'BEGIN {{OFS=FS}} {{ $1 = "\\""$1"\\""; print }}\' "{file_path}" > "{file_path}.tmp" && mv "{file_path}.tmp" "{file_path}"'])


### make a coverageBed file from the merged files to be used to normalise CRAC reads
### this file will be computed to RPK
# for file in glob.glob(bam_folder + '/WT*ALL_merged.bam'):
#     os.system('coverageBed -s -sorted -counts -a {} -b {} > {}/{}.txt'.format(reference_gff, file, CoverageBed_folder, os.path.basename(file)[:-4]))


# for file in glob.glob(CoverageBed_folder + '/WT*ALL_merged.txt'):
#     df = pd.read_csv(file, sep = '\t', header = None)
#     RPK_list = []
#     for index, row in df.iterrows():
#         RPK_list.append(int(row[9] * 1000 / (abs(row[3] - row[4]))))
    
#     df[9] = RPK_list
#     df.to_csv(file + '_RPK.txt', sep = '\t', header = None, index = False)

    
# for file in glob.glob(CoverageBed_folder + '/WT*ALL_merged.txt_RPK.txt'):
#     os.system('cut -f9-10 {} > {}_cut.txt'.format(file, file[:-4]))
#     os.system('rm {}'.format(file))


# # ### Only keep the name in the first column
# for file in glob.glob(CoverageBed_folder + '/WT*ALL_merged.txt_RPK_cut.txt'):
#     names = []
#     df = pd.read_csv(file, sep = '\t', header = None)
#     for item in df[0].tolist():
#         alist = item.split(';')
#         x = 0
#         for a in alist:
#             if a.find('gene_name=') == 0 or a.find('name=') == 0:
#                 names.append(a[10:])
#                 x += 1
#             else:
#                 pass
#         if x == 0:
#             print(item)
#             sys.exit()

#     df[0] = names

#     df.to_csv(file, sep = '\t', header = None, index = False)


# ### DESeq2 can't handle gene names containing space characters (' ')
# ### Unless it has quotation marks around
# for file_path in glob.glob(CoverageBed_folder + "/WT*ALL_merged.txt_RPK_cut.txt"):
#     subprocess.run(["bash", "-c", f'awk -F"\t" \'BEGIN {{OFS=FS}} {{ $1 = "\\""$1"\\""; print }}\' "{file_path}" > "{file_path}.tmp" && mv "{file_path}.tmp" "{file_path}"'])
