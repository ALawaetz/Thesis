import os
import glob
import pandas
import sys



motif_file = sys.argv[1]
motif = sys.argv[2]
output_file = f'{motif_file[:-4]}_{motif}.gtf'
select_motifs_script = '/Users/andreas/Bacillus/Edinburgh/pyCRAC/pycrac/pyCRAC/scripts/pySelectMotifsFromGTF.py'
bash_script = '/Users/andreas/Bacillus/Edinburgh/pyCRAC/pycrac/pyCRAC/Modify_extracted_motif_GTF_file.sh'


os.system(f'python {select_motifs_script} --gtf {motif_file} -m {motif} -o {output_file}')

os.system(f'bash {bash_script} {output_file}')