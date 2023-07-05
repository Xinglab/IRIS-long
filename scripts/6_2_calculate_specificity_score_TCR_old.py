import os,re,sys
from collections import defaultdict
import numpy as np

TCR_bound_peptide_inf_name = sys.argv[1]
protein_inf = sys.argv[2]
window_size = int(sys.argv[3])
required_HLA_str = sys.argv[4]
CPM_inf_name = sys.argv[5]
tumor_sample_count = int(sys.argv[6])
outf_dir = sys.argv[7]

######## find transcripts encoding peptides bound by given HLA complex ###########
peptide_dict = defaultdict()
peptide2HLA_dict = defaultdict()
with open(TCR_bound_peptide_inf_name, 'r') as inf_2:
	for index, line in enumerate(inf_2):
		arr = line.strip().split('\t')
		if index == 0: continue
		HLA_list = re.sub("\*","",arr[1]).split(';')
		for each_HLA in HLA_list:
			if re.findall(each_HLA, required_HLA_str):
				peptide2HLA_dict[arr[0]] = arr[1]
				peptide_dict[arr[0]] = []
				break
	
trans_ID = ''
with open(protein_inf, 'r') as inf_3:
	for line in inf_3:
		if line.startswith('>'):
			trans_ID = line.split('|')[2].split('_')[0]
		else:
			line = line.strip()
			for i in range(len(line)-window_size+1):
				this_peptide = line[i:i+window_size].upper()
				if this_peptide in peptide2HLA_dict:
					peptide_dict[this_peptide].append(trans_ID)
print("Transcripts contain sliding peptides are found.")

outf1_name = outf_dir+'/6_2_peptide_matrix_WindowSize_'+str(window_size)+'_TCR.txt'
outf1 = open(outf1_name,'w')
outf1.write('Peptide\tBound_by_HLA\tDerived_transcripts\n')
for each_peptide in sorted(peptide_dict.keys()):
	outf1.write(each_peptide+'\t'+peptide2HLA_dict[each_peptide]+'\t'+';'.join(peptide_dict[each_peptide])+'\n')
outf1.close()


################ Calculate tumor-specificity score#################

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return raw_ID.split('.')[0]+'-PAR-Y'
	else:
		return raw_ID.split('.')[0]

sample_list = []
exp_dict = defaultdict(lambda: [])
with open(CPM_inf_name, 'r') as CPM_inf:
	for index,line in enumerate(CPM_inf):
		arr = line.strip().split('\t')
		if index == 0:
			sample_list = arr[3:len(arr)]
			continue
		else:
			trans_ID = get_right_ID(arr[0])
			gene_ID = get_right_ID(arr[2])
			exp_dict[trans_ID] = np.array(list(map(float, arr[3:len(arr)])))

###### output sample list for each group #######
outf_sample = open(outf_dir+"/6_2_Sliding_sample_list.txt", 'w')
outf_sample.write('Tumor\t'+';'.join(sample_list[0:tumor_sample_count])+'\n')
outf_sample.write('Normal_tissue\t'+';'.join(sample_list[tumor_sample_count:len(sample_list)])+'\n')
outf_sample.close()
################################################
#peptide2exp = defaultdict(lambda: defaultdict(lambda: 0))
full_matrix = np.zeros(shape=(1,len(sample_list)))
full_matrix_log2 = np.matrix(['Peptide']+sample_list)

file_total_line = float(os.popen("wc -l %s" % outf1_name).readlines()[0].split(' ')[0])
with open(outf1_name, 'r') as inf_1:
	for index,line in enumerate(inf_1):
		arr = line.strip().split('\t')
		if index == 0: continue
		if index % 10000 == 0: 
			print(" %s of specificity score calculation is finished." % (str(round(index*100/file_total_line, 2))+'%'))
		this_peptide = arr[0]
		only_exp_list = np.zeros(len(sample_list))
		for each_trans in arr[2].split(';'):
			if len(exp_dict[each_trans]) == 0: continue
			only_exp_list += exp_dict[each_trans]
		full_matrix = np.vstack((full_matrix, only_exp_list))
		only_log2_list = [this_peptide] + list(map(str, np.log2(np.array(only_exp_list)+1)))
		full_matrix_log2 = np.vstack((full_matrix_log2, only_log2_list))

pseudo_count = 0.1
print ('tumor_sample', tumor_sample_count)
print ('normal sample', len(sample_list) - tumor_sample_count)
full_matrix_2 = full_matrix.transpose()
tumor_mean_row = np.mean(full_matrix_2[0:tumor_sample_count,1:], axis=0)
tissue_mean_row = np.mean(full_matrix_2[tumor_sample_count:len(sample_list),1:], axis=0)
FC_mean_row = np.log2((np.array(tumor_mean_row)+pseudo_count) / (np.array(tissue_mean_row)+pseudo_count))

tumor_median_row = np.median(full_matrix_2[0:tumor_sample_count,1:], axis=0)
tissue_median_row = np.median(full_matrix_2[tumor_sample_count:len(sample_list),1:], axis=0)
FC_median_row = np.log2((np.array(tumor_median_row)+pseudo_count) / (np.array(tissue_median_row)+pseudo_count))

full_matrix_log2 = full_matrix_log2.transpose()
#full_matrix_log2_rep = full_matrix_log2
full_matrix_log2 = np.vstack((full_matrix_log2, ['Log2FC_mean']+list(map(str, FC_mean_row))))
full_matrix_log2 = np.vstack((full_matrix_log2, ['Log2FC_median']+list(map(str, FC_median_row))))
full_matrix_log2 = full_matrix_log2.transpose()
outf_name = outf_dir+"/6_2_specificity_score_matrix_WindowSize_%s_TCR.txt" % (str(window_size))
with open(outf_name,'w') as outf:
	for line in full_matrix_log2:
		np.savetxt(outf, line, delimiter='\t', fmt='%s')

