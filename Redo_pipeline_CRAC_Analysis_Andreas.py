import os
import sys
import glob
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns
import math
# import BinCollector_analysis_Andreas
import shutil
import gff3


### Complete pipeline for analysis from novo files downwards

root_folder = sys.argv[1]
HTF_genes = sys.argv[2].split(',')

### Save version of script in root folder
script_location = os.path.realpath(__file__)
shutil.copy(script_location, f'{root_folder}/{os.path.basename(script_location)}')

### gtf file that has been adjusted for RNAseq and CRAC analysis by
### running the pyRAP_Github_publication.gff3 file through the scripts:
### 'exclude_all_tRNAs_rRNAs_and_associated.py', 'Prepare_pyRAP_to_CRAC_7June2023', and 'omit_overlapping_annotations.py'
### use the script Prepare_pyRAP_to_CRAC_7June2023.py to make CRAC friendly GTF file from pyRAP GFF file
### use the script exclude_all_tRNAs_rRNAs_and_associated.py to make the gtf file
### To remove tRNAs and rRNAs run this code after the above commands
# import gff3
# file = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files/B_subtilis_pyRAP_Github_publication.gff3_no_rRNAs_tRNAs_or_associated.gff3_CRAC_modified_no_overlaps.gtf'
# df = pd.read_csv(file, sep = '\t', names = gff3.header)
# df = df[(df['source'] != 'tRNA') & (df['source'] != 'rRNA') & (df['source'] != 'rRNA_or_tRNA_adjacent_transcript') & (df['source'] != 'Novel_rRNA_or_tRNA_adjacent')]
# output_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files'
# df.to_csv(os.path.join(output_folder, 'B_subtilis_CRAC_RNAseq_analysis.gtf'), sep = '\t', header = None, index = False)

# file = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files/B_subtilis_pyRAP_Github_publication.gff3_no_rRNAs_tRNAs_or_associated.gff3_CRAC_modified._no_overlaps.gff3'
# df = pd.read_csv(file, sep = '\t', names = gff3.header)
# df = df[(df['source'] != 'tRNA') & (df['source'] != 'rRNA') & (df['source'] != 'rRNA_or_tRNA_adjacent_transcript') & (df['source'] != 'Novel_rRNA_or_tRNA_adjacent')]
# output_folder = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files'
# df.to_csv(os.path.join(output_folder, 'B_subtilis_CRAC_RNAseq_analysis.gff3'), sep = '\t', header = None, index = False)


gtf_file = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Annotation_files/B_subtilis_CRAC_RNAseq_analysis.gtf'
chrom_file = '/Users/andreas/Bacillus/Bioinformatics/Txt_and_word_files/B_subtilis_chromsize.txt'
tab_file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_30thMay2024/B_subtilis_NC_000964.3.fasta.tab'
coverageBed_file = '/Users/andreas/Bacillus/Bioinformatics/RNAseq_jag_khpa_kre/X204SC22071444-Z01-F001/Analysis/Coverage_bed/WT_3_FKRN220334190-1A_H3J2FDSX5_L1_ALL_merged.txt_RPK_cut.txt'
GTF2Bedgraph_script = '/Users/andreas/Bacillus/Edinburgh/pyCRAC/pycrac/pyCRAC/scripts/pyGTF2bedGraph.py'
figure_folder = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures'

os.makedirs(root_folder + '/DPW', exist_ok = True)
os.makedirs(root_folder + '/pyReadCounters_RPKM_normalised_to_RNAseq', exist_ok = True)
os.makedirs(root_folder + '/1d_pyReadCounters_RPKM_normalised_to_RNAseq', exist_ok = True)

### Make all the folders
for HTF_gene in HTF_genes:
    HTF_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis'
    os.makedirs(HTF_folder + '/1_pyReadCounters', exist_ok = True)
    os.makedirs(HTF_folder + '/1b_pyReadCounters_wig_deletions', exist_ok = True)
    os.makedirs(HTF_folder + '/1c_pyGTF2Bedgraph', exist_ok = True)
    os.makedirs(HTF_folder + '/2_pyCalculateFDRs', exist_ok = True)
    os.makedirs(HTF_folder + '/3_Significant_reads', exist_ok = True)
    os.makedirs(HTF_folder + '/4_pyClusterReads', exist_ok = True)
    os.makedirs(HTF_folder + '/4_pyClusterReads_delsOnly', exist_ok = True)
    os.makedirs(HTF_folder + '/4_pyClusterReads_delsOnly_Targets_only', exist_ok = True)
    os.makedirs(HTF_folder + '/5_Significant_reads_overlapping_Clusters', exist_ok = True)
    os.makedirs(HTF_folder + '/5b_Significant_reads_overlapping_Clusters_wig_deletions', exist_ok = True)
    os.makedirs(HTF_folder + '/6_pyReadCounters_signficant_reads_overlapping_clusters', exist_ok = True)
    os.makedirs(HTF_folder + '/7_TPMPW', exist_ok = True)

 
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

################################################### pyReadCounters on novo files ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyReadCounters {HTF_gene}')
#     input_folder = f'{root_folder}/novo_files_merged'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1_pyReadCounters'
#     for file in glob.glob(input_folder + f'/{HTF_gene}*.novo'):
#         if os.path.basename(file).find('WT') == -1:
#             os.system('python pyReadCounters.py -f {} --gtf {} --sense --blocks --mutations=delsonly -o {}/{}'.format(file, gtf_file, output_folder, os.path.basename(file)))



############################################# Calculate DPW ##############################################################
## Deletions pr WT
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Calculate DPW for {HTF_gene}')
#     # input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5b_Significant_reads_overlapping_Clusters_wig_deletions'
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1b_pyReadCounters_wig_deletions'
#     output_folder = input_folder
#     for file in glob.glob(input_folder + '/*.wig'):
#         os.system(f'python DPW.py {file} {gtf_file} {HTF_gene} {output_folder}')
#     ### merge the scoreboard files
#     df_merge = pd.DataFrame()
#     for file in glob.glob(input_folder + '/*_scoreboard.txt'):
#         df = pd.read_csv(file, sep = '\t', header = 0)
#         df_merge = pd.concat([df_merge, df])
#     df_merge = df_merge.sort_values(['DPW'], ascending = False)
#     df_merge.to_csv(os.path.join(f'{root_folder}/DPW', f'{HTF_gene}_fwd_and_rev_scoreboard.txt'), sep = '\t', index = False)




############################################ pyGTF2Bedgraph on pyReadCounters files ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyGTF2Bedgraph on pyReadCounters files {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1_pyReadCounters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1c_pyGTF2Bedgraph'
#     for file in glob.glob(input_folder + f'/{HTF_gene}*.gtf'):
#         if os.path.basename(file).find('WT') == -1:
#             output_file = f'{output_folder}/{os.path.basename(file)[:-4]}'
#             os.system(f'python {GTF2Bedgraph_script} --gtf {file} -o {output_file} -c {chrom_file}')

# ### adjust bedgraph files so there are no overlapping lines (compatability with bedGraphToBigWig)
# for HTF_gene in HTF_genes:
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1c_pyGTF2Bedgraph'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1c_pyGTF2Bedgraph'
#     for file in glob.glob(input_folder + f'/{HTF_gene}*.bedgraph'):
#         if os.path.basename(file).find('WT') == -1:
#             df = pd.read_csv(file, sep = '\t', header = None)
#             df[2] = [i - 1 for i in df[2].tolist()]
#             df.to_csv(file, sep = '\t', header = None, index = False)

# ### Now convert to bigwig         
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'converting bedgraph files to bigwig {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1c_pyGTF2Bedgraph'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1c_pyGTF2Bedgraph'
#     for file in glob.glob(input_folder + f'/{HTF_gene}*.bedgraph'):
#         if os.path.basename(file).find('WT') == -1:
#             output_file = f'{output_folder}/{os.path.basename(file)[:-9]}'
#             os.system(f'bedGraphToBigWig {file} {chrom_file} {output_file}.bw')


# ############################################# Make wig deletions files on all reads ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Make wig deletion files for {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1_pyReadCounters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1b_pyReadCounters_wig_deletions'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system(f'python deletions_wig_files_pipeline.py {file} {output_folder} {HTF_gene} {chrom_file}')



# # ####################################################### pyCalculateFDRs ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyCalculateFDRs {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1_pyReadCounters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/2_pyCalculateFDRs'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system('python pyCalculateFDRs.py -f {} --gtf {} -o {}/{} -c {} -m 0.001'.format(file, gtf_file, output_folder, os.path.basename(file), chrom_file))


# # # ######################################### Bedtools intersect pyReadCounters and pyCalculateFDRs #############################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Bedtools intersect {HTF_gene} for finding FDR signficant reads')
#     pyReadCounters = f'{root_folder}/{HTF_gene}_CRAC_analysis/1_pyReadCounters'
#     pyCalculateFDRs = f'{root_folder}/{HTF_gene}_CRAC_analysis/2_pyCalculateFDRs'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/3_Significant_reads'
#     for a_file, b_file in zip(sorted(glob.glob(pyReadCounters + '/*.gtf')), sorted(glob.glob(pyCalculateFDRs + '/*.gtf'))):
#         print('Intersecting {}___{}'.format(os.path.basename(a_file), os.path.basename(b_file)))
#         os.system('bedtools intersect -s -u -a {} -b {} > {}/{}'.format(a_file, b_file, output_folder, os.path.basename(a_file)))


# # # # ####################################################### pyClusterReads ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyClusterReads {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/3_Significant_reads'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system('python pyClusterReads.py -f {} --gtf {} --cic 10 --co 10 --ch 50 -o {}/{} --mutsfreq 10'.format(file, gtf_file, output_folder, os.path.basename(file)))



# # # # # # # ############################################### Remove clusters without X% deletions #############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Removing clusters without deletions {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system(f'bash clusters_muts_only.sh {file} {output_folder}')


# # # # # # # ############################## Merge overlapping annotions in cluster files #############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Merging annotations in cluster file {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system(f'python merge_annotations_Andreas.py {file}')


# # # # # ######################################### Bedtools intersect Significant Reads and Clusters #############################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Bedtools intersect {HTF_gene} for finding signficant reads in clusters')
#     significant_reads = f'{root_folder}/{HTF_gene}_CRAC_analysis/3_Significant_reads'
#     clusters = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5_Significant_reads_overlapping_Clusters'
#     for a_file, b_file in zip(sorted(glob.glob(significant_reads + '/*.gtf')), sorted(glob.glob(clusters + '/*.gtf'))):
#         print('Intersecting {}___{}'.format(os.path.basename(a_file), os.path.basename(b_file)))
#         os.system('bedtools intersect -s -u -a {} -b {} > {}/{}'.format(a_file, b_file, output_folder, os.path.basename(a_file)))


# # ############################################ Make wig deletions files on significant reads ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Make wig deletion files for {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5_Significant_reads_overlapping_Clusters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5b_Significant_reads_overlapping_Clusters_wig_deletions'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system(f'python deletions_wig_files_pipeline.py {file} {output_folder} {HTF_gene} {chrom_file}')


############################################# pyReadCounters on significant clusters ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyReadCounters {HTF_gene} on significant clusters')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5_Significant_reads_overlapping_Clusters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/6_pyReadCounters_signficant_reads_overlapping_clusters'
#     for file in glob.glob(input_folder + f'/{HTF_gene}*.gtf'):
#         if os.path.basename(file).find('WT') == -1:
#             os.system('python pyReadCounters.py -f {} --gtf {} --file_type gtf --sense -o {}/{}'.format(file, gtf_file, output_folder, os.path.basename(file)))


# ########################################## Calculate RPKM and normalise to RNAseq ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Calculate RPKM and normalise to RNAseq for {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/6_pyReadCounters_signficant_reads_overlapping_clusters'
#     output_folder = root_folder + '/1d_pyReadCounters_RPKM_normalised_to_RNAseq'
#     cluster_height = 50
#     for file in glob.glob(input_folder + f'/{HTF_gene}*hittable_cDNAs.txt'):
#         os.system(f'python RPKM_Andreas.py {file} {gtf_file} {output_folder} {coverageBed_file} {cluster_height}')



# ########################################## Order genes by strongest cross-linking ####################################################################
# ### Check the script Order_top_CRAC_hits.py in the pycrac folder to see if all file names are correct
# ### set the minimum ratio for something to be a hit
# the_ratio = 10
# the_cutoff = 0.01
# os.system(f"python Order_top_CRAC_hits.py {the_ratio} {the_cutoff}")

######################################## Prepare list for GO enrichment analysis ####################################################################
# Check the script Prepare_CRAC_for_GO_analysis in the pycrac folder to see if all file names are correct
os.system(f"python Prepare_CRAC_for_GO_analysis.py")


# ################################# Remove clusters that are not on the list of "most cross-linked genes" ####################################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Remove non-target clusters {HTF_gene}')
#     cluster_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly'
#     ratio_folder = root_folder + '/ratios'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly_Targets_only'
#     for cluster_file, ratio_file in zip(glob.glob(cluster_folder + f'/{HTF_gene}*.gtf'), glob.glob(ratio_folder + f'/{HTF_gene}*.txt')):
#         print(cluster_file + '.........' + ratio_file)
#         os.system(f'python Keep_Target_Clusters_Only.py {cluster_file} {ratio_file} {output_folder}')


######################################################## pyMotif ####################################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyMotif {HTF_gene}')
#     features = ['All', 'CDS', '5UTR', '3UTR']
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly'
#     for file in glob.glob(input_folder + '/*final_merge'):
#         for feature in features:
#             os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_dels_{feature}', exist_ok = True)
#             output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_dels_{feature}'
#             if feature == 'All':
#                 os.system('python pyMotif.py -f {} --file_type gtf --gtf {} --tab {} -o {}/{} --k_min 7 --k_max 7'.format(file, gtf_file, tab_file, output_folder, os.path.basename(file)))
#             else:
#                 os.system('python pyMotif.py -f {} --file_type gtf --gtf {} --tab {} -o {}/{} -a {} --k_min 7 --k_max 7'.format(file, gtf_file, tab_file, output_folder, os.path.basename(file), feature))
       


########################################### Calculate RPKM and normalise to RNAseq ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Calculate RPKM and normalised to RNAseq for {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/6_pyReadCounters_signficant_reads_overlapping_clusters'
#     output_folder = root_folder + '/pyReadCounters_RPKM_normalised_to_RNAseq'
#     for file in glob.glob(input_folder + f'/{HTF_gene}*.txt'):
#         os.system(f'python RPKM_Andreas.py {file} {gtf_file} {output_folder} {coverageBed_file}')




############################################# Calculate DPW ##############################################################
## Deletions pr WT
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Calculate DPW for {HTF_gene}')
#     # input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5b_Significant_reads_overlapping_Clusters_wig_deletions'
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/1b_pyReadCounters_wig_deletions'
#     output_folder = input_folder
#     for file in glob.glob(input_folder + '/*.wig'):
#         os.system(f'python DPW.py {file} {gtf_file} {HTF_gene} {output_folder}')
#     ### merge the scoreboard files
#     df_merge = pd.DataFrame()
#     for file in glob.glob(input_folder + '/*_scoreboard.txt'):
#         df = pd.read_csv(file, sep = '\t', header = 0)
#         df_merge = pd.concat([df_merge, df])
#     df_merge = df_merge.sort_values(['DPW'], ascending = False)
#     df_merge.to_csv(os.path.join(f'{root_folder}/DPW', f'{HTF_gene}_fwd_and_rev_scoreboard.txt'), sep = '\t', index = False)


# ############################################ Plot distribution of features ##############################################################
# print('###########################')
# print('Plotting distribution of features')
# input_folder = f'{root_folder}/ratios'
# top_genes = 200
# os.system(f'python feature_distribution.py {input_folder} {gtf_file} {top_genes} {figure_folder}')


# ############################################ pyBinCollector ##############################################################
# features = ['All', 'CDS', '5UTR', '3UTR']
# # features = ['CDS']
# os.chdir('/Users/andreas/Bacillus/Edinburgh/pyCRAC/pycrac/pyCRAC')
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyBinCollector {HTF_gene}')
#     reads_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/5_Significant_reads_overlapping_Clusters/{HTF_gene}_merged.novo_count_output_cDNAs.gtf'
#     cluster_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_delsOnly/{HTF_gene}_merged.novo_count_output_cDNAs.gtf'
#     ratios_file = f'{root_folder}/ratios/{HTF_gene}_ratios.txt'
#     top_genes = 800
#     for feature in features:
#         os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyBinCollector_{feature}', exist_ok = True)
#         output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyBinCollector_{feature}'
#         os.system(f'python pyBinCollector_pipeline_21stJune2024.py {reads_file} {cluster_file} {gtf_file} {root_folder} {output_folder} {feature} {ratios_file} {top_genes}')


# ############################################ Add features to ratios files ##############################################################
# os.makedirs(f'{root_folder}/ratios/features/', exist_ok = True)
# output_folder = f'{root_folder}/ratios/features/'
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Add features to ratios files {HTF_gene}')
#     file = f'{root_folder}/ratios/{HTF_gene}_ratios.txt'
#     df = pd.read_csv(file, sep = '\t', header = 0)
#     gtf = pd.read_csv(gtf_file, sep = '\t', names = gff3.header)
#     gtf_names = [i.split(' ')[1] for i in gtf['attributes'].tolist()]
#     gtf['name'] = gtf_names
#     df_features = []
#     for item in df['Gene'].tolist():
#         df_features.append(gtf['source'][gtf['name'] == item].iloc[0])
#     df['feature'] = df_features
#     df.to_csv(os.path.join(output_folder, os.path.basename(file)), sep = '\t', header = True, index = False)


# # # ########################################### pyReadCounters on significant reads in clusters ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyReadCounters {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/5_Significant_reads_overlapping_Clusters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/6_pyReadCounters_signficant_reads_overlapping_clusters'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         os.system('python pyReadCounters.py -f {} --file_type gtf --gtf {} --sense --blocks --mutations=delsonly -o {}/{}'.format(file, gtf_file, output_folder, os.path.basename(file)))


# #################################################### Calculate TPMPW ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Calculate transcripts pr million pr WT RNAseq (TPMPW) {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/6_pyReadCounters_signficant_reads_overlapping_clusters'
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/7_TPMPW'
    
#     if len(glob.glob(input_folder + '/*.txt')) > 1:
#         print('Write new code, this only works when replicates have been merged')
#         sys.exit()
#     else:
#         hittable_file = glob.glob(input_folder + '/*.txt')[0]
    
#     if HTF_gene == 'SpoVG':
#         os.system(f'python TPMPW.py {hittable_file} {output_folder} SpoVG')
#     else:
#         os.system(f'python TPMPW.py {hittable_file} {output_folder} not_SpoVG')

    ### This block of code can be used if a non-pyRAP annotation file is used
    # if HTF_gene == 'SpoVG':
    #     os.system(f'python TPMPW_SpoVG.py {hittable_file} {output_folder}')
    # else:
    #     os.system(f'python TPMPW.py {hittable_file} {output_folder}')


################################################ Find top 200 clusters ##############################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Finding top 200 clusters {HTF_gene}')
#     cluster_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads'
#     top200 = pd.read_csv(glob.glob(f'{root_folder}/{HTF_gene}_CRAC_analysis/7_TPMPW' + '/*.txt')[0], sep = '\t', header = 0)['Gene'].tolist()[:200]
#     file = f'{cluster_folder}/{HTF_gene}_merged.novo_count_output_cDNAs.gtf'
#     # Define a custom function to filter rows starting with '#'
#     def skip_comment_rows(row):
#         return row.startswith('#')
#     df = pd.read_csv(file, sep = '\t', header = None, skiprows=4)
#     names = []
#     for attribute in df[8].tolist():
#         attribute_list = attribute.split('; ')
#         for a in attribute_list:
#             asplit = a.split(' ')
#             if asplit[0] == 'gene_name':
#                 names.append(asplit[1].replace('"', ''))
#                 break
#             else:
#                 pass
    
#     index_list = []
#     for i in range(len(names)):
#         gene_list = names[i].split(',')
#         for gene in gene_list:
#             if gene in top200:
#                 index_list.append(i)
#                 break
#             else:
#                 pass
#     df = df[df.index.isin(index_list)]
#     df.to_csv(f'{file[:-4]}_top200.gtf', sep = '\t', index = False, header = None)

                

######################################################## pyMotif ####################################################################
# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'pyMotif {HTF_gene}')
#     features = ['All', 'CDS', '5UTR', '3UTR']
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads_dels'
#     for file in glob.glob(input_folder + '/*.gtf'):
#         for feature in features:
#             os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_dels_{feature}', exist_ok = True)
#             output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_dels_{feature}'
#             if feature == 'All':
#                 os.system('python pyMotif.py -f {} --file_type gtf --gtf {} --tab {} -n 100 -o {}/{} --k_min 6 --k_max 8'.format(file, gtf_file, tab_file, output_folder, os.path.basename(file)))
#             else:
#                 os.system('python pyMotif.py -f {} --file_type gtf --gtf {} --tab {} -n 100 -o {}/{} -a {} --k_min 6 --k_max 8'.format(file, gtf_file, tab_file, output_folder, os.path.basename(file), feature))
        

############################################# assign gene categories to TPMPW ####################################################################
### To assign categories to CRAC hits
# gene_categories = '/Users/andreas/Bacillus/Bioinformatics/CRAC_16Jan2024/gene_categories_Andreas.txt'
# gene_categories = pd.read_csv(gene_categories, sep = '\t', header = 0, dtype = 'str')
# gene_categories = gene_categories.fillna('NA')

# thegroups = []
# for column in gene_categories.columns.tolist():
#     if column != 'locus' and column != 'title' and column != 'synonyms' and column != 'outlinks->bsu' and column != 'outlinks->bsupath':
#         thegroups.append(column)

# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Assign gene categories to TPMPW {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/DPW'
#     os.makedirs(f'{input_folder}_categories', exist_ok = True)
#     output_folder = f'{input_folder}_categories'

#     for file in glob.glob(input_folder + '/*.txt'):
#         df = pd.read_csv(file, sep = '\t', header = 0)
#         df = df[df.index < 100]
#         df_split = pd.DataFrame()
#         names = []
#         TPMPW = []
#         for index, row in df.iterrows():
#             thegene = row['Gene'].replace('BSU_misc_RNA_', 'BSU@misc@RNA@')
#             thegene = thegene.replace('BSU_', 'BSU@')
#             gene_list = thegene.split('_')
#             gene_list = [i.replace('BSU@misc@RNA@', 'BSU_misc_RNA_') for i in gene_list]
#             gene_list = [i.replace('BSU@', 'BSU_') for i in gene_list]
#             for gene in gene_list:
#                 names.append(gene)
#                 TPMPW.append(row['DPW'])
#         df_split['Gene'] = names
#         df_split['DPW'] = TPMPW
#         del df
#         df = df_split

#         new_df = pd.DataFrame()

#         unassigned = []
#         categories = []
#         bcategories = []
#         for thegene in df['Gene'].tolist():
#             alist = []
#             blist = []
#             if thegene.find('Novel_transcript') == 0:
#                 alist.append('NA')
#                 blist.append('NA')
#             else:
#                 thegene = thegene.replace('BSU_misc_RNA_', 'BSU@misc@RNA@')
#                 thegene = thegene.replace('BSU_', 'BSU@')
#                 gene_list = thegene.split('_')
#                 gene_list = [i.split('.')[0] for i in gene_list]
#                 gene_list = [i.replace('BSU@misc@RNA@', 'BSU_misc_RNA_') for i in gene_list]
#                 gene_list = [i.replace('BSU@', 'BSU_') for i in gene_list]

#                 for gene in gene_list:
#                     if gene == 'miscRNA':
#                         alist.append('NA')
#                         blist.append('NA')
#                     else:
#                         for group in thegroups:
#                             try:
#                                 if gene_categories[group][gene_categories['title'] == gene].iloc[0] != 'NA':
#                                     alist.append('{}: {}'.format(group, gene_categories[group][gene_categories['title'] == gene].iloc[0]))
#                                     blist.append(group)
#                                 else:
#                                     alist.append('NA')
#                                     blist.append('NA')
#                             except:
#                                 unassigned.append(gene)

#             astring = ', '.join(list(set([i for i in alist if i != 'NA'])))
#             bstring = ', '.join(list(set([i for i in blist if i != 'NA'])))
#             categories.append(astring)
#             bcategories.append(bstring)
#         df['Categories'] = categories
#         df.to_csv(os.path.join(output_folder, f'{os.path.basename(file)[:-4]}_categories.txt'), sep = '\t', index = False)
#         df['Categories'] = bcategories
#         df.to_csv(os.path.join(output_folder, f'{os.path.basename(file)[:-4]}_categories_groups.txt'), sep = '\t', index = False)
#         df = pd.DataFrame()
#         df['unassigned'] = unassigned
#         df.to_csv(os.path.join(output_folder, f'{os.path.basename(file)[:-4]}_unassigned.txt'), sep = '\t', index = False)





# ################################################# fisher exact test ####################################################################

# ### for fisher's exact test we need a two-by-two table
# ###                          =Category                    !=Category
# ### InGroup                 InGroup_Cat                   InGroup_NonCat
# ### NonGroup                NonGroup_Cat                  NonGroup_NonCat       


# gene_categories_revised = '/Users/andreas/Bacillus/Bioinformatics/CRAC_16Jan2024/gene_categories_Andreas_revised.txt'
# gene_categories_revised = pd.read_csv(gene_categories_revised, sep = '\t', header = 0, dtype = 'str')
# gene_categories_revised = gene_categories_revised.fillna('NA')

# convert_categories = '/Users/andreas/Bacillus/Bioinformatics/CRAC_16Jan2024/categories-2023-05-22.csv'
# convert_categories = pd.read_csv(convert_categories, sep = ',', header = 0)

# file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_16Jan2024/gene_categories_Andreas.txt'
# converter = pd.read_csv(file, sep = '\t', usecols = ['locus', 'title'])

# ### How many genes to analyse
# top_number_of_genes = 500

# for HTF_gene in HTF_genes:
#     the_source = 'CRAC'
#     print('###########################')
#     print(f'fisher exact test {HTF_gene}')
#     output_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/DPW_categories'

#     file = f'{output_folder}/{HTF_gene}_fwd_and_rev_scoreboard_categories_groups.txt'
#     df_CRAC = pd.read_csv(file, sep = '\t', header = 0)
#     df_CRAC = df_CRAC.dropna().reset_index(drop = True)

#     ### Add Locus to dataframes
#     locus_list = []
#     for name in df_CRAC['Gene'].tolist():
#         name = name.replace('@', ' ')
#         name = name.split('.')[0]
#         if name.find('Novel_transcript') == 0:
#             locus_list.append('Novel_transcript')
#         else:
#             try:
#                 locus_list.append(converter['locus'][converter['title'] == name].iloc[0])
#             except:
#                 print('problem here A')
#                 print(name)
#                 sys.exit()

#     df_CRAC['Locus'] = locus_list

#     ### For more accurate statistics we will omit Novel_transcripts as they are not part of the
#     ### SubtiWiki gene_categories_revised dataset
#     df_CRAC = df_CRAC[df_CRAC['Locus'] != 'Novel_transcript']

#     ### We will also remove duplicate rows so e.g. dnaA.1 and dnaA.2 are not both counted
#     df_CRAC = df_CRAC.drop_duplicates(['Locus'])
#     df_CRAC = df_CRAC.reset_index(drop = True)
#     df_CRAC = df_CRAC[df_CRAC.index < top_number_of_genes]


#     thecategories = []
#     for category in df_CRAC['Categories'].tolist():
#         alist = category.split(', ')
#         alist = [i for i in alist if i != 'NoCategories' and i != 'NA']
#         thecategories += alist
#     thecategories = list(set(thecategories))

#     sig_df = pd.DataFrame()
#     sig_category = []
#     sig_p_value = []
#     sig_odds_ratio = []

#     ### Subset gene_categories_revised so only contains genes that are not in df_CRAC
#     merged = gene_categories_revised.merge(df_CRAC, on='Locus', how='left', indicator=True)
#     if the_source == 'CRAC':
#         gene_categories_revised = merged[merged['_merge'] == 'left_only'].drop(columns=['_merge', 'Gene', 'Categories_y', 'DPW'])
#         gene_categories_revised = gene_categories_revised.rename(columns={'Categories_x': 'Categories'})
#     else:
#         gene_categories_revised = merged[merged['_merge'] == 'left_only'].drop(columns=['_merge', 'Gene', 'Categories_y'])
#         gene_categories_revised = gene_categories_revised.rename(columns={'Categories_x': 'Categories'})


#     for category in thecategories:
#         InGroup_Cat = 0
#         InGroup_NonCat = 0
#         NonGroup_Cat = 0
#         NonGroup_NonCat = 0

#         for item, locus in zip(df_CRAC['Categories'].tolist(), df_CRAC['Locus'].tolist()):
#             categories_row = item.split(', ')
#             if category in categories_row:
#                 InGroup_Cat += 1
#             elif category not in categories_row:
#                 InGroup_NonCat += 1
#             else:
#                 print('something wrong, review code')
        

#         for item, locus in zip(gene_categories_revised['Categories'].tolist(), gene_categories_revised['Locus'].tolist()):
#             categories_row = item.split(', ')
#             if category in categories_row:
#                 NonGroup_Cat += 1
#             elif category not in categories_row:
#                 NonGroup_NonCat += 1
#             else:
#                 print('something wrong, review code')

        
#         ### We'll use alternative=Greater
#         ### Because we are only interested in gene categories that are enriched not the oppposite
#         table = np.array([[InGroup_Cat, InGroup_NonCat], [NonGroup_Cat, NonGroup_NonCat]])
#         res = fisher_exact(table, alternative = 'greater')
#         odds = res[0]
#         pvalue = res[1]
#         if pvalue <= 0.05:
#             sig_category.append(category)
#             sig_p_value.append(pvalue)
#             sig_odds_ratio.append(odds)


#     new_cat = []
#     for item in sig_category:
#         new_cat.append(convert_categories['category'][convert_categories['id'] == item].iloc[0])

#     sig_df['Category'] = new_cat
#     sig_df['P_value'] = sig_p_value
#     sig_df['Odds_ratio'] = sig_odds_ratio
#     sig_df = sig_df.sort_values(['Odds_ratio', 'P_value'], ascending = False)


#     sig_df.to_csv(os.path.join(output_folder, f'{HTF_gene}_significant_categories.csv'), sep = '\t', index = False)





# ################################################# Distribution of features ####################################################################

# print('###########################')
# print(f'Distribution of features plot')
# figure_folder = '/Users/andreas/Bacillus/Bioinformatics/CRAC_25Sep23/figures'
# pyRAP = pd.read_csv(gtf_file, sep = '\t', header = None)
# names = []
# for item in pyRAP[8].tolist():
#     alist = item.split(';')
#     names.append(alist[0][10:])
# pyRAP['name'] = names
# df_stack = pd.DataFrame()
# df_stack['feature'] = list(set(pyRAP[1].tolist()))
# for HTF_gene in HTF_genes:
#     df_HTF_list = []
#     htf_basefolder = f'{root_folder}/{HTF_gene}_CRAC_analysis'
#     TPMPW_folder = f'{htf_basefolder}/7_TPMPW'
#     df = pd.read_csv(glob.glob(TPMPW_folder + '/*.txt')[0], sep = '\t', header = 0)
#     df = df[df.index < 200]
    
#     top100genes = df['Gene'].tolist()
#     if HTF_gene == 'SpoVG':
#         print(top100genes)
#     roots = []
#     for gene in top100genes:
#         gene_feature = pyRAP[1][pyRAP['name'] == gene].iloc[0]
#         if (gene.split('.')[0], gene_feature) in roots:
#             pass
#         else:
#             roots.append((gene.split('.')[0], gene_feature))
            
#     lengths = []
#     for feature in df_stack['feature'].tolist():
#         lengths.append(len([i for i in roots if i[1] == feature]))
#     lengths = [i * 100 / sum(lengths) for i in lengths]
#     df_stack[HTF_gene] = lengths


# numeric_columns = df_stack.select_dtypes(include=['int', 'float']).columns
# df_stack['Sum'] = df_stack[numeric_columns].sum(axis=1)
# df_stack = df_stack[df_stack['Sum'] != 0]
# df_stack = df_stack.sort_values(['Sum'], ascending = False)
# df_stack = df_stack.drop('Sum', axis = 1)
# df_stack = df_stack.set_index('feature')



# ### set figure parameters
# SMALL_SIZE = 8
# MEDIUM_SIZE = 10
# BIGGER_SIZE = 12
# LARGE_SIZE = 14


# plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=LARGE_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# plt.rc('font', family='sans-serif')  # Specify a different font family

# df_stack = df_stack.transpose()

# # Plot the stacked bar chart
# fig, ax = plt.subplots()

# ax = df_stack.plot.bar(stacked=True)

#     # Set the x-axis tick positions
# ax.set_xticks(range(len(df_stack.index)))
# # Set the x-axis tick labels and rotate them by 90 degrees
# ax.set_xticklabels(df_stack.index, rotation = 45)

# # Adding labels and title
# ax.set_xlabel('Genes')
# ax.set_ylabel('%')
# ax.set_title('Distribution of features')
# ax.legend()

# # Reverse the order of the legend handles and labels
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc='upper left')

# # plt.show()

# ### save plot
# plt.savefig(os.path.join(figure_folder, 'Distribution_of_features.pdf'), bbox_inches='tight')
# df_stack.to_csv(os.path.join(figure_folder, 'top200genes_features.csv'), sep = '\t')


################################################# Heatmeaps ####################################################################

# print('###########################')
# print(f'Heatmap plot')
# df_merge = pd.DataFrame()
# df_merge['Gene'] = []
# for HTF_gene in HTF_genes:
#     TPMPW = f'{root_folder}/{HTF_gene}_CRAC_analysis/7_TPMPW'
#     df = pd.read_csv(glob.glob(TPMPW + '/*.txt')[0], sep = '\t', header = 0, names = ['Gene', HTF_gene])#, comment = '#', usecols = [0, 1], names = ['Gene', 'Hits'])
#     df = df[df.index < 200]
#     df[HTF_gene] = [i * 10 for i in df[HTF_gene].tolist()] ### to avoid negative values
#     df[HTF_gene] = df[HTF_gene].apply(lambda x: math.log10(x) if x != 0 else x)
#     df_merge = pd.merge(df_merge, df, on='Gene', how='outer')

# df_merge = df_merge.fillna(0)
# df_merge = df_merge.set_index('Gene')


# # Plot the heatmap
# plt.figure(figsize=(8, 6))
# heatmap = sns.heatmap(df_merge, cmap='YlOrRd', cbar_kws={"shrink": 0.3})

# # Get the y-axis tick positions
# if len(df_merge) > 50:
#     plt.yticks([])
# else:
#     yticks = range(len(df_merge.index))
#     yticks_adjusted = [tick + 0.5 for tick in yticks]

#     # Set the y-axis tick positions and labels
#     plt.yticks(yticks_adjusted, df_merge.index)

# plt.xticks(rotation=45)

# # Set the axis labels)
# plt.ylabel('Genes')

# # Set the title
# plt.title('Gene hits', fontsize = 16)

# # Set the legend title
# cbar = heatmap.collections[0].colorbar
# cbar.ax.set_title('Log10')

# # Display the plot
# plt.savefig(os.path.join(figure_folder, 'Heatmap_TPMPW.pdf'), bbox_inches='tight')
# plt.close()




######################################## plot only deletions in clusters BinCollector ####################################################################
### To assign categories to CRAC hits

# for HTF_gene in HTF_genes:
#     print('###########################')
#     print(f'Extracting dels coordinates from clusters {HTF_gene}')
#     input_folder = f'{root_folder}/{HTF_gene}_CRAC_analysis/4_pyClusterReads'
#     file = f'{input_folder}/{HTF_gene}_merged.novo_count_output_cDNAs_top200.gtf'

#     features = ['All'] + ['CDS', '5UTR', '3UTR', 'sRNA_independent', 'transcript', 'sRNA_5UTR', 'intergenic_UTR', 'sRNA_3UTR', 'putative_sRNA', 'S_feature', 'sRNA_inter_intragenic', 'Novel_transcript', 'pseudogene']
#     small_features = ['5UTR', '3UTR', 'sRNA_5UTR', 'intergenic_UTR', 'sRNA_3UTR', 'sRNA_inter_intragenic']
#     for feature in features:
#         os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/10_pyBinCollector_dels', exist_ok = True)
#         output_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/10_pyBinCollector_dels/{os.path.basename(file)[:-4]}_{feature}.txt'
#         os.system(f'python pyBinCollector.py -f {file} --gtf {gtf_file} -a {feature} -n 50 -o {output_file} --deletions --normalize --outputall')

# for HTF_gene in HTF_genes:
#     input_file = f'/Users/andreas/Bacillus/Bioinformatics/CRAC_analysis_17June/{HTF_gene}_CRAC_analysis/10_pyBinCollector_dels/{HTF_gene}_merged.novo_count_output_cDNAs_top200_CDS.txt'
#     plot_title = f'{HTF_gene}_CDS'
#     figure_file_name = plot_title

#     BinCollector_analysis_Andreas.BinCollector_plot(input_file, plot_title, figure_folder, figure_file_name)


#     del_df = pd.DataFrame()
#     df = pd.read_csv(file, sep = '\t', header = None, skiprows = 4)
#     start = []
#     percent = []
#     strand_list = []
#     new_attributes = []
#     for attribute, strand in zip(df[8].tolist(), df[6].tolist()):
#         attribute_list = attribute.split('# ')
#         deletions = attribute_list[-1].replace(';', '')
#         deletions = deletions.split(',')
#         for deletion in deletions:
#             deletion = deletion.split('D')
#             start.append(int(deletion[0]))
#             percent.append(deletion[1])
#             strand_list.append(strand)
#             new_attributes.append(attribute)
    
#     del_df[0] = ['NC_000964.3'] * len(start)
#     del_df[1] = ['cluster'] * len(start)
#     del_df[2] = ['interval'] * len(start)
#     del_df[3] = start
#     del_df[4] = [i + 1 for i in start]
#     del_df[5] = percent
#     del_df[6] = strand_list
#     del_df[7] = ['.'] * len(start)
#     del_df[7] = new_attributes
#     output_file = f'{input_folder}/{os.path.basename(file)[:-4]}_dels_coordinates.gtf'
#     del_df.to_csv(os.path.join(input_folder, f'{os.path.basename(file)[:-4]}_dels_coordinates.gtf'), sep = '\t', header = None, index = False)

#     features = ['All'] + ['CDS', '5UTR', '3UTR', 'sRNA_independent', 'transcript', 'sRNA_5UTR', 'intergenic_UTR', 'sRNA_3UTR', 'putative_sRNA', 'S_feature', 'sRNA_inter_intragenic', 'Novel_transcript', 'pseudogene']
#     small_features = ['5UTR', '3UTR', 'sRNA_5UTR', 'intergenic_UTR', 'sRNA_3UTR', 'sRNA_inter_intragenic']
#     for feature in features:
#         os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/10_pyBinCollector_dels', exist_ok = True)
#         output_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/10_pyBinCollector_dels/{os.path.basename(file)[:-4]}_{feature}.txt'
#         if feature in small_features:
#             os.system(f'python pyBinCollector.py -f {output_file} --gtf {gtf_file} -a {feature} -n 20 -o {output_file}')
#         else:
#             os.system(f'python pyBinCollector.py -f {output_file} --gtf {gtf_file} -a {feature} -n 50 -o {output_file}')