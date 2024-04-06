import os,re,sys
from collections import defaultdict
import numpy as np

TCR_bound_peptide_inf_name = sys.argv[1]
protein_inf = sys.argv[2]
window_size = int(sys.argv[3])
required_HLA_str = sys.argv[4]
CPM_inf_name = sys.argv[5]
tumor_sample_count = int(sys.argv[6])
cds_inf = sys.argv[7]
outf_dir = sys.argv[8]


def get_geome_loc(chrom, strand, peptide_loc, cds):  ### ('chr1', '+','1-9','56629545#56629592;56631205#56631690')
	seq_coords_list = []
	pep2nt_start = (int(peptide_loc.split('-')[0]) - 1) * 3
	pep2nt_end = int(peptide_loc.split('-')[1]) * 3 - 1   #### 1-9 in pep pos => 0-26 in nt pos
	[find_start_tag_global, find_end_tag_global] = [0, 0]
	current_nt_pos = 0
	if strand == "+":
		for each_cds_seg in cds.split(";"):
			[find_start_tag, find_end_tag] = [0, 0]
			[each_cds_start, each_cds_end] = list(map(int, each_cds_seg.split("#")))
			each_cds_len = each_cds_end - each_cds_start + 1
			if current_nt_pos <= pep2nt_start <= current_nt_pos + each_cds_len:   ## pep starts in this exon
				this_start_pos = each_cds_start + (pep2nt_start - current_nt_pos)
				[find_start_tag, find_start_tag_global] = [1, 1]
			if current_nt_pos <= pep2nt_end <= current_nt_pos + each_cds_len:   ## pep ends in this exon
				this_end_pos = each_cds_start + (pep2nt_end - current_nt_pos)
				[find_end_tag, find_start_tag_global] = [1, 1]
			current_nt_pos += each_cds_len

			if (find_start_tag == 0) and (find_end_tag == 0):
				if find_start_tag_global == 1:
					seq_coords_list.append(f"{each_cds_start}-{each_cds_end}")
			elif (find_start_tag == 1) and (find_end_tag == 0):
				seq_coords_list.append(f"{this_start_pos}-{each_cds_end}")
			elif (find_start_tag == 0) and (find_end_tag == 1):
				seq_coords_list.append(f"{each_cds_start}-{this_end_pos}")
				break
			elif (find_start_tag == 1) and (find_end_tag == 1):
				seq_coords_list.append(f"{this_start_pos}-{this_end_pos}")
				break
		seq_coords = '#'.join(seq_coords_list)

	elif strand == "-":   ### 59035010#59035126;59036471#59036573;59036690#59036745
		for each_cds_seg in cds.split(";")[::-1]:
			[find_start_tag, find_end_tag] = [0, 0]
			[each_cds_start, each_cds_end] = list(map(int, each_cds_seg.split("#")))
			each_cds_len = each_cds_end - each_cds_start + 1
			#### since it's negative strand, current_nt_pos actually means the number of nt from 3' to 5' direction
			if current_nt_pos <= pep2nt_start <= current_nt_pos + each_cds_len:   ## pep starts in this exon
				this_start_pos = each_cds_end - (pep2nt_start - current_nt_pos)
				[find_start_tag, find_start_tag_global] = [1, 1]
			if current_nt_pos <= pep2nt_end <= current_nt_pos + each_cds_len:   ## pep ends in this exon
				this_end_pos = each_cds_end - (pep2nt_end - current_nt_pos)
				[find_end_tag, find_start_tag_global] = [1, 1]
			current_nt_pos += each_cds_len

			if (find_start_tag == 0) and (find_end_tag == 0):
				if find_start_tag_global == 1:
					seq_coords_list.append(f"{each_cds_start}-{each_cds_end}")
			elif (find_start_tag == 1) and (find_end_tag == 0):
				seq_coords_list.append(f"{each_cds_start}-{this_start_pos}")
			elif (find_start_tag == 0) and (find_end_tag == 1):
				seq_coords_list.append(f"{this_end_pos}-{each_cds_end}")
				break
			elif (find_start_tag == 1) and (find_end_tag == 1):
				seq_coords_list.append(f"{this_end_pos}-{this_start_pos}")
				break
		seq_coords = '#'.join(seq_coords_list[::-1])

	genome_coords = "_".join([chrom, strand, seq_coords])
	return genome_coords


######## make nucleotide position dictionary ########
nucleotide_dict = defaultdict()
with open(cds_inf, "r") as cds_file:
	for index,line in enumerate(cds_file):
		if index == 0: continue
		arr = line.strip().split("\t")
		transcriptID = arr[0]
		nucleotide_dict[transcriptID] = [arr[3], arr[4], arr[5]]  ## [chromsome, strand, cds]


######## find transcripts encoding peptides bound by given HLA complex ###########
peptide_dict = defaultdict()
peptide2HLA_dict = defaultdict()
with open(TCR_bound_peptide_inf_name, 'r') as inf_2:
	for index, line in enumerate(inf_2):
		arr = line.strip().split('\t')
		if index == 0: continue
		HLA_list = re.sub(r"\*","",arr[1]).split(';')
		for each_HLA in HLA_list:
			if re.findall(each_HLA, required_HLA_str):
				peptide2HLA_dict[arr[0]] = arr[1]
				peptide_dict[arr[0]] = {"trans_ID" : [], "peptide_pos": [], "genomic_pos": []}
				break
	
trans_ID = ''
with open(protein_inf, 'r') as inf_3:
	for line in inf_3:
		if line.startswith('>'):
			trans_ID = line.split(' ')[0].split("|")[-1].split('_')[0]
		else:
			line = line.strip()
			[chrom, strand, cds] = nucleotide_dict[trans_ID]
			for i in range(len(line)-window_size+1):
				this_peptide = line[i:i+window_size].upper()
				if this_peptide in peptide2HLA_dict:
					peptide_loc = f"{i+1}-{i+window_size}"
					genome_loc = get_geome_loc(chrom, strand, peptide_loc, cds)
					# update the dictionary
					peptide_dict[this_peptide]["trans_ID"].append(trans_ID)
					peptide_dict[this_peptide]["peptide_pos"].append(peptide_loc)
					peptide_dict[this_peptide]["genomic_pos"].append(genome_loc)
					
print("Transcripts contain sliding peptides are found.")

outf1_name = outf_dir+'/6_2_peptide_matrix_WindowSize_'+str(window_size)+'_TCR.txt'
outf1 = open(outf1_name,'w')
outf1.write('Peptide\tBound_by_HLA\tDerived_transcripts\tPeptide_position\tGenomic_position\n')
for each_peptide in sorted(peptide_dict.keys()):
	outf1.write(each_peptide+'\t'+peptide2HLA_dict[each_peptide]+'\t'+';'.join(peptide_dict[each_peptide]["trans_ID"])+'\t'+';'.join(peptide_dict[each_peptide]["peptide_pos"])+'\t'+';'.join(peptide_dict[each_peptide]["genomic_pos"])+'\n')
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
full_matrix_log2 = np.matrix(['Peptide']+sample_list)

file_total_line = float(os.popen("wc -l %s" % outf1_name).readlines()[0].split(' ')[0])
with open(outf1_name, 'r') as inf_1:
	for index,line in enumerate(inf_1):
		arr = line.strip().split('\t')
		if index == 0: continue
		if index % 1000 == 0: 
			print(" %s of specificity score calculation is finished." % (str(round(index*100/file_total_line, 2))+'%'))
		this_peptide = arr[0]
		only_exp_list = np.zeros(len(sample_list))
		for each_trans in arr[2].split(';'):
			if len(exp_dict[each_trans]) == 0: continue
			only_exp_list += exp_dict[each_trans]
		if index == 1:
			full_matrix = np.matrix(only_exp_list).reshape(1,len(only_exp_list))
		else:
			full_matrix = np.vstack((full_matrix, only_exp_list))
		only_log2_list = [this_peptide] + list(map(str, np.log2(np.array(only_exp_list)+1)))
		full_matrix_log2 = np.vstack((full_matrix_log2, only_log2_list))

pseudo_count = 0.1
print ('tumor_sample', tumor_sample_count)
print ('normal sample', len(sample_list) - tumor_sample_count)

tumor_mean_col = np.mean(full_matrix[0:,0:tumor_sample_count], axis=1)
tissue_mean_col = np.mean(full_matrix[0:,tumor_sample_count:len(sample_list)], axis=1)
FC_mean = np.log2((np.array(tumor_mean_col)+pseudo_count) / (np.array(tissue_mean_col)+pseudo_count))
FC_mean_col = np.matrix(FC_mean).reshape(len(FC_mean),1)

tumor_median_col = np.median(full_matrix[0:,0:tumor_sample_count], axis=1)
tissue_median_col = np.median(full_matrix[0:,tumor_sample_count:len(sample_list)], axis=1)
FC_median = np.log2((np.array(tumor_median_col)+pseudo_count) / (np.array(tissue_median_col)+pseudo_count))
FC_median_col = np.matrix(FC_median).reshape(len(FC_median),1)

full_matrix_log2 = np.hstack((full_matrix_log2, np.vstack((['Log2FC_mean'], FC_mean_col))))
full_matrix_log2 = np.hstack((full_matrix_log2, np.vstack((['Log2FC_median'], FC_median_col))))

outf_name = outf_dir+"/6_2_specificity_score_matrix_WindowSize_%s_TCR.txt" % (str(window_size))
with open(outf_name,'w') as outf:
	for line in full_matrix_log2:
		np.savetxt(outf, line, delimiter='\t', fmt='%s')
