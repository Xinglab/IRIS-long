import os,re,sys
from collections import defaultdict

window_size = int(sys.argv[4])
target_trans = sys.argv[1]
protein_inf_name = sys.argv[2]
out_dir = sys.argv[3]
outf_name = out_dir+'/7_1_peptide_matrix_WindowSize_'+str(window_size)+'_'+target_trans+'.txt'

trans_ID = ''
peptide_dict = defaultdict(lambda: [])
peptide_pos_dict = defaultdict(lambda: [])
with open(protein_inf_name, 'r') as inf_2:
	for line in inf_2:
		if line.startswith('>'):
			trans_ID = line.split('|')[2].split('_')[0]
		else:
			line = line.strip()
			if trans_ID == target_trans: 
				for i in range(len(line)-window_size+1):
					this_peptide = line[i:i+window_size].upper()
					peptide_dict[this_peptide].append(trans_ID)
					peptide_pos_dict[this_peptide].append(str(i+1))
				break
print("Sliding peptides from target protein are built.")

with open(protein_inf_name, 'r') as inf_3:
	for line in inf_3:
		if line.startswith('>'):
			trans_ID = line.split('|')[2].split('_')[0]
		else:
			line = line.strip()
			for i in range(len(line)-window_size+1):
				this_peptide = line[i:i+window_size].upper()
				if this_peptide in peptide_dict:
					if trans_ID not in peptide_dict[this_peptide]:
						peptide_dict[this_peptide].append(trans_ID)
print("Transcripts contain sliding peptides are found.")

outf = open(outf_name,'w')
outf.write('Peptide\tPos_in_protein\tDerived_transcripts\n')
for each_peptide in sorted(peptide_dict.keys()):
	outf.write(each_peptide+'\t'+';'.join(peptide_pos_dict[each_peptide])+'\t'+';'.join(peptide_dict[each_peptide])+'\n')
outf.close()

