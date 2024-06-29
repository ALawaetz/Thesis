import os
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
import BinCollector_analysis_Andreas



root_folder = sys.argv[1]
HTF_gene = sys.argv[2]
the_motif = sys.argv[3]
the_feature = sys.argv[4]

### simple txt file with list of genes
### The code searches for the motif in these genes if they are in the cluster file and makes the plot on those genes
### example genes_file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_Edinburgh_Jag_KhpA_Kre_SpoVG/genes_antimicrobial_toxin_phage_TA.txt'
genes_file = sys.argv[5]

select_motifs_script = '/Users/andreas/Bacillus/Edinburgh/pyCRAC/pycrac/pyCRAC/scripts/pySelectMotifsFromGTF.py'
bash_script_cure_MotifSelect = '/Users/andreas/Bacillus/Edinburgh/pyCRAC/pycrac/pyCRAC/Modify_extracted_motif_GTF_file.sh'

figure_folder = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/figures'


# ####################################################################################################################################
### pySelectMotifFromGTF 

def motif_bin_plot(HTF_gene, motif, feature, genes_file):
    print('###########################')
    print(f'pySelectMotifFromGTF {HTF_gene}')
    motif_file = glob.glob(f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_dels_{feature}' + '/*.gtf')[0]
    thelength = len(motif)
    os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_{feature}_MotifSelect', exist_ok = True)
    SelectMotif_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/8_pyMotif_{feature}_MotifSelect/{os.path.basename(motif_file)[:-4]}_{motif}.gtf'
    Reads_file = glob.glob(f'{root_folder}/{HTF_gene}_CRAC_analysis/5_Significant_reads_overlapping_Clusters' + '/*.gtf')[0]
    os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/9_pyBinCollector', exist_ok = True)
    BinCollector_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/9_pyBinCollector/{os.path.basename(SelectMotif_file)}'
    thesource = f'{feature}'
    thefeature = f'{feature}'

    os.system(f'python {select_motifs_script} --gtf {motif_file} -m {motif} -o {SelectMotif_file}')
    os.system(f'bash {bash_script_cure_MotifSelect} {SelectMotif_file} {thesource} {thefeature}')
    

    # # ### pyBinCollector with Select Motifs
    # # ### HTF_gene
    os.makedirs(f'{root_folder}/{HTF_gene}_CRAC_analysis/9_pyBinCollector', exist_ok = True)
    if thesource == 'All':
        os.system(f'python pyBinCollector.py -f {Reads_file} --gtf {SelectMotif_file} -n {100 + thelength} -r 50 --max_length {100 + thelength + 2} --deletions --outputall --normalize -o {BinCollector_file}')
    else:
        os.system(f'python pyBinCollector.py -f {Reads_file} --gtf {SelectMotif_file} -a {thesource} -n {100 + thelength} -r 50 --max_length {100 + thelength + 2} --deletions --outputall --normalize -o {BinCollector_file}')
    sys.exit()
    # ### Plot BinCollector results
    TPMPW_file = f'{root_folder}/{HTF_gene}_CRAC_analysis/7_TPMPW/{HTF_gene}_merged.novo_count_output_cDNAs.gtf_hittable_cDNAs.txt'
    BinCollector_analysis_Andreas.BinCollector_plot(BinCollector_file, f'{HTF_gene} {thesource} {motif}', figure_folder, f'{HTF_gene}_{motif}_BinCollector_plot', genes_file)

motif_bin_plot(HTF_gene, the_motif, the_feature, genes_file)





# def individual_bin_plot(Gene, file, save, motif):
#     df = pd.read_csv(file, sep = '\t', header = 0)
#     # Get the row with index name 'X'
#     row = df.loc['murAA.2']
#     column_names = row.index
#     y_values = row.values
#     plt.bar(column_names, y_values)
#     plt.xlabel('Bp')
#     plt.ylabel('Read_intensity')
#     plt.title(f'{Gene} {motif}')
#     plt.xticks(rotation=90)
#     plt.xticks(range(0, len(column_names), 10))
#     if save == 'save':
#         plt.savefig(os.path.join(figure_folder, f'{Gene}_{motif}.pdf'), bbox_inches='tight')
#     plt.show()
#     plt.close()

# file = '/Users/andreas/Bacillus/Bioinformatics/CRAC_jag_khpA_kre_19thJune2024/KhpA_CRAC_analysis/9_pyBinCollector/KhpA_merged.novo_count_output_cDNAs.gtffinal_merge_CDS_top_k-mers_in_features_YTGCCGS.gtf'
# Gene = 'murAA.2'
# individual_bin_plot(Gene, file, 'save', 'YTGCCGS')






