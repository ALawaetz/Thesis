import os
import glob
import shutil
import pandas as pd
import subprocess
import sys

bam_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/bam_files'
bam_folder_stranded = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/bam_files_stranded'
reference_gff = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_spoVG_LB_OD3/Analysis/B_subtilis_pyRAP_no_header.gff3'
CoverageBed_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed'
bedgraph_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/bedgraph_files'
norm_bedgraph_all_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/bedgraph_files_norm_all'
norm_bedgraph_jag_khpA_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/bedgraph_files_norm_jag_khpA'
norm_bedgraph_kre_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/bedgraph_files_norm_kre'
chrom_size = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/B_subtilis_chromsize.txt'

for folder in [bam_folder_stranded, CoverageBed_folder, bedgraph_folder, norm_bedgraph_all_folder, norm_bedgraph_jag_khpA_folder, norm_bedgraph_kre_folder]:
    if not os.path.exists(folder):
        os.makedirs(folder)
    

### Merge bam with the same genotype
### We can then make bedgraph files of these files to view in Jbrowse
### Instead of having three tracks with WT and three tracks with mutant 1 etc. we
### will have just one track pr. genotype
### We can normalise coverage by multiplying each count in the bedgraph files
### with the estimateSizeFactors() function in DESeq2 in R.
# bam_dict = {}
# for file in glob.glob(bam_folder + '/*.bam'):
#     if 'khpA_1' in os.path.basename(file):
#         pass
#     elif os.path.basename(file)[:2] in list(set(bam_dict.keys())):
#         bam_dict[os.path.basename(file)[:2]] += [file]
#     else:
#         bam_dict[os.path.basename(file)[:2]] = [file]


# for key, item in bam_dict.items():
#     print('merging files in {}: {}'.format(key, item))
#     print('\n')
#     if len(item) == 3:
#         os.system('samtools merge {}/{}_ALL_merged.bam {} {} {}'.format(bam_folder, os.path.basename(item[0])[:-4], item[0], item[1], item[2]))
#     elif len(item) == 2:
#         os.system('samtools merge {}/{}_ALL_merged.bam {} {}'.format(bam_folder, os.path.basename(item[0])[:-4], item[0], item[1]))
#     else:
#         print('More or less than three replicates. Review code.')
#         break


### sort bam files in forward and reverse reads
# for file in glob.glob(bam_folder + '/*.bam'):
#     print(f'making fwd and rev strand files for {os.path.basename(file)}')
#     os.system('samtools view -bS -F 16 {} > {}/{}_fwd_strand.bam'.format(file, bam_folder_stranded, os.path.basename(file[:-4])))
#     os.system('samtools view -bS -f 16 {} > {}/{}_rev_strand.bam'.format(file, bam_folder_stranded, os.path.basename(file[:-4])))


### Make bedgraph coverage files and convert to BigWig
# os.system('cut -f1,4-5,9,6,7 {} > {}.bed'.format(reference_gff, reference_gff[:-5]))
# for file in glob.glob(bam_folder_stranded + '/*strand.bam'):
#     print(f'bedtools intersect on {os.path.basename(file)}')
#     os.system('bedtools intersect -a {} -b {}.bed > {}_intersected.bam'.format(file, reference_gff[:-5], file[:-4]))


# for file in glob.glob(bam_folder_stranded + '/*intersected.bam'):
#     print(f'making bedgraph file with {os.path.basename(file)}')
#     os.system('bedtools genomecov -trackline -bg -ibam {} > {}.bedgraph'.format(file, file[:-4]))


### Before converting the bedgraph files to BigWig
### We will first normalise according to "read depth"
### This is done best by using estimateSizeFactors(dds) in DESeq2
### Reference: median ratio method Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 2010, 11:106. http://dx.doi.org/10.1186/gb-2010-11-10-r106
### This is also explained here https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
## sort bam files
# for file in glob.glob(bam_folder + '/*ALL_merged.bam'):
#     print(f'sorting {os.path.basename(file)}')
#     os.system('samtools sort {} -T {} -o {}'.format(file, os.path.basename(file), file))

# ### make coverageBed files to input DESeq2
# for file in glob.glob(bam_folder + '/*ALL_merged.bam'):
#     os.system('coverageBed -s -sorted -counts -a {} -b {} > {}/{}.txt'.format(reference_gff, file, CoverageBed_folder, os.path.basename(file)[:-4]))

# for file in glob.glob(CoverageBed_folder + '/*ALL_merged.txt'):
#     os.system('cut -f9-10 {} > {}_cut.txt'.format(file, file[:-4]))
#     os.system('rm {}'.format(file))


# ### Only keep the name in the first column
# for file in glob.glob(CoverageBed_folder + '/*ALL_merged_cut.txt'):
#     print(os.path.basename(file))
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


### Move bedgraph files to bedgraph folder
# for file in glob.glob(bam_folder_stranded + '/*.bedgraph'):
#     if 'khpA_1' in os.path.basename(file):
#         pass
#     else:
#         os.system(f'mv {file} {bedgraph_folder}/{os.path.basename(file)}')


### divide value column in bedgraph files with size factor
### size factor are calculated in R using DESeq2 and the estimateSizeFactors(dds)
### Size factors vary depending on what samples are included. If strains are ['jag', 'khpA', 'kre', 'WT'] or only ['jag', 'WT']
### Therefore compute size factors depending on the purpose and the comparison which is to be made
def norm_bedgraph(the_size_factors, the_strains, the_bedgraph_folder, the_norm_folder, search, exclude):
    size_factor_dict = {}
    for sf, s in zip(the_size_factors, the_strains):
        size_factor_dict[s] = sf
    
    print('size_factor_dict:')
    print(size_factor_dict)

    # Get the list of all files matching the criteria in search and exclude
    filtered_files = glob.glob(the_bedgraph_folder + f'/*{search}*.bedgraph')
    for ex in exclude:
        filtered_files = [f for f in filtered_files if not (ex in os.path.basename(f))]

    for file in sorted(filtered_files):
        for key, item in size_factor_dict.items():
            if os.path.basename(file).find(key) == 0:
                the_size_factor = item
                break
            else:
                pass
        
        print(f'normalising {os.path.basename(file)} using sizefactor {the_size_factor}')

        text = open(file, 'r').readlines()
        newtext = text[0]
        for item in text[1:]:
            alist = item.rstrip('\n').split('\t')
            astring = ''
            astring += '{}\t'.format(alist[0])
            astring += '{}\t'.format(alist[1])
            astring += '{}\t'.format(alist[2])
            astring += '{}\n'.format(int(int(alist[3]) / the_size_factor))
            newtext += astring
        norm_bedgraph_file = the_norm_folder + f'/{os.path.basename(file[:-9])}_normalised.bedgraph'
        with open(norm_bedgraph_file, 'w') as f:
            f.writelines(newtext)
            f.close()

    for file in glob.glob(the_norm_folder + '/*.bedgraph'):
        os.system('bedGraphToBigWig {} {} {}.bw'.format(file, chrom_size, file))


# ### Size factors when all strains are merged according to genotype and included in the same comparison
# # size_factors = [1.1722519, 0.6825686, 1.2902603, 0.9727421]
# # strains = ['jag', 'khpA', 'kre', 'WT']
# # norm_folder = norm_bedgraph_all_folder
# # norm_bedgraph(size_factors, strains, bedgraph_folder, norm_folder, 'ALL_merged', [''])


# ### Size factors when all strains are included in the same comparison and NOT merged according to genotype  
# size_factors = [1.1740485, 0.8514995, 1.1632301, 0.9036818, 0.9518536, 1.0758802, 1.1812295, 1.2313488, 0.8320279, 0.8244540, 0.9847106]
# strains = ['jag_1', 'jag_2', 'jag_3', 'khpA_2', 'khpA_3', 'kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3']
# norm_folder = norm_bedgraph_all_folder
# norm_bedgraph(size_factors, strains, bedgraph_folder, norm_folder, '', ['ALL_merged'])


# ### Size factors when strains (jag, khpA, WT) are merged according to genotype and included in the same comparison
# size_factors = [1.2785469, 0.7390437, 1.0563001]
# strains = ['jag', 'khpA', 'WT']
# norm_folder = norm_bedgraph_jag_khpA_folder
# norm_bedgraph(size_factors, strains, bedgraph_folder, norm_folder, 'ALL_merged', ['kre'])

# ### Size factors when strains (jag, khpA, WT) are included in the same comparison and NOT merged according to genotype and 
# size_factors = [1.2446468, 0.8991003, 1.2208945, 0.9543690, 1.0013825, 0.8780895, 0.8665175, 1.0310542]
# strains = ['jag_1', 'jag_2', 'jag_3', 'khpA_2', 'khpA_3', 'WT_1', 'WT_2', 'WT_3']
# norm_folder = norm_bedgraph_jag_khpA_folder
# norm_bedgraph(size_factors, strains, bedgraph_folder, norm_folder, '', ['ALL_merged', 'kre'])


# ### Size factors when strains (jag, khpA, WT) are merged according to genotype and included in the same comparison
size_factors = [1.151580, 0.868372]
strains = ['kre', 'WT']
norm_folder = norm_bedgraph_kre_folder
norm_bedgraph(size_factors, strains, bedgraph_folder, norm_folder, 'ALL_merged', ['jag', 'khpA'])


# ### Size factors when strains (jag, khpA, WT) are included in the same comparison and NOT merged according to genotype and 
size_factors = [1.0627624, 1.1764206, 1.2204356, 0.8263015, 0.8199746, 0.9786852]
strains = ['kre_1', 'kre_2', 'kre_3', 'WT_1', 'WT_2', 'WT_3']
norm_folder = norm_bedgraph_kre_folder
norm_bedgraph(size_factors, strains, bedgraph_folder, norm_folder, '', ['ALL_merged', 'jag', 'khpA'])
