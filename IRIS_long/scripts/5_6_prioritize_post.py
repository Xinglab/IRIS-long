import os,re,sys
from collections import defaultdict
import numpy as np

tumor_count = int(sys.argv[1])
cutoff_score = float(sys.argv[2]) # 1
isoform_abundance_inf_name = sys.argv[3] 
pc_fasta_inf_name = sys.argv[4] # ../3_PARCB_PC.fasta
cds_inf = sys.argv[5] # ./4.4_match_ID.txt
window_size = int(sys.argv[6]) # default: 9
outf_dir = sys.argv[7] # ./

pseudo_count = 0.1
tissue_count = 0

######## make nucleotide position dictionary ########
nucleotide_dict = defaultdict()
with open(cds_inf, "r") as cds_file:
	for index,line in enumerate(cds_file):
		if index == 0: continue
		arr = line.strip().split("\t")
		transcriptID = arr[0]
		nucleotide_dict[transcriptID] = [arr[3], arr[4], arr[5]]  ## [chromsome, strand, cds]

#### load candidate transcripts ####
required_trans_list = []
with open("%s/5_5_Summarized_CAR_T_prioritized_targets_temp.txt" % outf_dir, "r") as cand_trans_inf:
	for index,line in enumerate(cand_trans_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		required_trans_list.append(arr[0])

#### obtain abundance of trans ####
tumor_cpm_dict = defaultdict(lambda: [])
tissue_cpm_dict = defaultdict(lambda: [])
with open(isoform_abundance_inf_name, 'r') as CPM_inf:
	for index,line in enumerate(CPM_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		trans_ID = arr[0]
		tumor_cpm_dict[trans_ID] = np.array(list(map(float, arr[3:3+tumor_count])))
		tissue_cpm_dict[trans_ID] = np.array(list(map(float, arr[3+tumor_count:len(arr)])))
		tissue_count = len(arr) - tumor_count - 3

#### match peptide to transcripts ####
pep2trans_dict = defaultdict(lambda: [])
trans2seq_dict = defaultdict()

trans_ID = ""
with open(pc_fasta_inf_name, "r") as protein_inf:
	for line in protein_inf:
		line = line.strip()
		if line.startswith(">"):
			trans_ID = line.split(' ')[0].split("|")[-1].split('_')[0]
		else:
			if trans_ID in required_trans_list:
				trans2seq_dict[trans_ID] = line
				for i in range(len(line)-window_size+1):
					this_peptide = line[i:i+window_size].upper()
					pep2trans_dict[this_peptide] = []

trans_ID = ""
with open(pc_fasta_inf_name, "r") as protein_inf2:
	for line in protein_inf2:
		line = line.strip()
		if line.startswith(">"):
			trans_ID = line.split(' ')[0].split("|")[-1].split('_')[0]
		else:
			for i in range(len(line)-window_size+1):
				this_peptide = line[i:i+window_size].upper()
				if this_peptide in pep2trans_dict:
					if trans_ID not in pep2trans_dict[this_peptide]:
						pep2trans_dict[this_peptide].append(trans_ID)

pep2score_dict = defaultdict()
for each_pep in pep2trans_dict:
	tumor_CPM_list = np.array([0.0]*tumor_count)
	tissue_CPM_list = np.array([0.0]*tissue_count)
	score = 0
	for each_trans in pep2trans_dict[each_pep]:
		if each_trans not in tumor_cpm_dict: continue
		tumor_CPM_list += tumor_cpm_dict[each_trans]
		tissue_CPM_list += tissue_cpm_dict[each_trans]
	score = round(np.log2((np.median(tumor_CPM_list)+pseudo_count)/(np.median(tissue_CPM_list)+pseudo_count)),2)
	if score >= cutoff_score:
		pep2score_dict[each_pep] = score

#### final output ####

#### function to concatenate strings ###
def concat(pos_list, score_list, protein_seq):
	conca_pos_list = []
	conca_score_list = []
	conca_pep_list = []
	pre_pos = pos_list[0]
	pre_score = score_list[0]
	for i in range(1,len(pos_list)):
		if float(score_list[i]) == float(pre_score):
			if int(pos_list[i].split('-')[0]) < int(pre_pos.split('-')[1]):
				pre_pos = pre_pos.split('-')[0]+'-'+pos_list[i].split('-')[1]
			else:
				conca_pos_list.append(pre_pos)
				conca_score_list.append(pre_score)
				pre_pos = pos_list[i]
				pre_score = score_list[i]
		else:
			conca_pos_list.append(pre_pos)
			conca_score_list.append(pre_score)
			pre_pos = pos_list[i]
			pre_score = score_list[i]
	conca_pos_list.append(pre_pos)
	conca_score_list.append(pre_score)
	for each_pos in conca_pos_list:
		this_pep = protein_seq[int(each_pos.split('-')[0])-1:int(each_pos.split('-')[1])]
		conca_pep_list.append(this_pep)
	
	#### combine a big range for regions with high specificity
	if len(conca_pos_list) == 1:
		final_conca_pos_list = conca_pos_list
	elif len(conca_pos_list) > 1:
		final_conca_pos_list = []
		pre_pos_final = conca_pos_list[0]
		for i in range(1,len(conca_pos_list)):
			if int(conca_pos_list[i].split('-')[0]) < int(pre_pos_final.split('-')[1]):
				pre_pos_final = pre_pos_final.split('-')[0]+'-'+conca_pos_list[i].split('-')[1]
			else:
				final_conca_pos_list.append(pre_pos_final)
				pre_pos_final = conca_pos_list[i]
		final_conca_pos_list.append(pre_pos_final)
	return [conca_pos_list, conca_score_list, conca_pep_list, final_conca_pos_list]

def get_overlap(specificity_pos_list, extra_pos_list, protein_seq):
	overlap_pos_list = []
	overlap_pep_list = []
	for each_region_1 in specificity_pos_list:
		[region_1_left, region_1_right] = list(map(int, each_region_1.split('-')))
		for each_region_2 in extra_pos_list:
			[region_2_left, region_2_right] = list(map(int, each_region_2.split('-')))
			if region_2_left >= region_2_right: continue
			if (region_2_right-region_1_left)*(region_2_left-region_1_right) < 0:
				overlap_left = max(region_1_left,region_2_left)
				overlap_right = min(region_1_right,region_2_right)
				overlap_pos_list.append(str(overlap_left)+'-'+str(overlap_right))
				overlap_pep_list.append(protein_seq[overlap_left-1:overlap_right])
	return(overlap_pos_list, overlap_pep_list)

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


outf = open("%s/5_5_Summarized_CAR_T_prioritized_targets_final.txt" % outf_dir, 'w')
with open("%s/5_5_Summarized_CAR_T_prioritized_targets_temp.txt" % outf_dir, "r") as cand_trans_inf_2:
	for index,line in enumerate(cand_trans_inf_2):
		line = line.strip()
		arr = line.split('\t')
		if index == 0: 
			outf.write(arr[0]+'\t'+arr[1]+'\t'+arr[11]+'\tTarget_peptide_pos\tTarget_peptide_seq\tTarget_peptide_in_genome_pos\t'+'\t'.join(arr[2:11])+'\t'+'\t'.join(arr[12:len(arr)])+'\tPeptide_with_high_tumor_specificity\tPos_pep_with_high_tumor_specificity\tTumor_specificity_score_for_each_regions\n')
			continue
		else:
			pep_list = []
			pos_list = []
			score_list = []
			protein_seq = trans2seq_dict[arr[0]]
			for i in range(len(protein_seq)-window_size+1):
				this_peptide = protein_seq[i:i+window_size].upper()
				if this_peptide in pep2score_dict:
					score_list.append(str(pep2score_dict[this_peptide]))
					pep_list.append(this_peptide)
					pos_list.append(str(i+1)+'-'+str(i+window_size))
			### candidates need to have at least one region exceed the cutoff of specificity score
			if len(score_list) > 0:
				if len(score_list) == 1:
					[concat_pos_list, concat_score_list, concat_pep_list, final_conca_pos_list] = [pos_list, score_list, pep_list, pos_list]
				else:
					[concat_pos_list, concat_score_list, concat_pep_list, final_conca_pos_list] = concat(pos_list, score_list, protein_seq)
				
				### candidates need to have overlap part between extracellular region and region with high specificity score
				if arr[-3] in ['Annotated_cell_surface:TM','Identified by both: consistent']:
					extra_string_pos = arr[-6].split(':')[1].split(',')
				else:
					extra_string_pos = ["1-"+str(len(protein_seq))]

				### print (final_conca_pos_list,extra_string_pos)
				[overlap_pos, overlap_pep] = get_overlap(final_conca_pos_list, extra_string_pos, protein_seq)

				### get corresponding genomic regions for final peptide sequence
				[chrom, strand, cds] = nucleotide_dict[arr[0]]
				genome_loc_list = []
				for peptide_loc in overlap_pos:
					genome_loc = get_geome_loc(chrom, strand, peptide_loc, cds)
					genome_loc_list.append(genome_loc)

				if len(overlap_pos) > 0:
					outf.write(arr[0]+'\t'+arr[1]+'\t'+arr[11]+'\t'+';'.join(overlap_pos)+'\t'+';'.join(overlap_pep)+'\t'+';'.join(genome_loc_list)+'\t'+'\t'.join(arr[2:11])+'\t'+'\t'.join(arr[12:len(arr)])+'\t'+';'.join(concat_pep_list)+'\t'+';'.join(concat_pos_list)+'\t'+';'.join(concat_score_list)+'\n')
outf.close()


