import os,re,sys
from collections import defaultdict, OrderedDict
import numpy as np

tumor_count = int(sys.argv[1])
tissue_CPM_cutoff = float(sys.argv[2])
tissue_percentage_cutoff = float(sys.argv[3])
isoform_cpm_inf_name = sys.argv[4]
isoform_prop_inf_name = sys.argv[5] 
annotated_isoform_contri_inf_name = sys.argv[6] 
tumor_vs_normal_inf_name = sys.argv[7]
gtf_inf_name = sys.argv[8]
outf_dir = sys.argv[9] # ./

other_isoform_prop_cutoff = 10
anno_percent_cutoff = 50
interested_isoform_prop_cutoff = 10

dir_path = os.path.dirname(os.path.realpath(__file__))
#### convert ID to gene name ####
def get_right_ID(raw_ID):
	if re.findall('_PAR_Y',raw_ID):
		new_ID = raw_ID.split('.')[0]+'-PAR-Y'
	else:
		new_ID = raw_ID.split('.')[0]
	return new_ID

ID2name_dict = defaultdict()
trans2gene_dict = defaultdict()

with open(gtf_inf_name, 'r') as gtf_inf:
	for line in gtf_inf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] == 'gene':
			gene_ID = get_right_ID(re.findall("gene_id \"(.+?)\";",arr[8])[0])
			gene_name = re.findall("gene_name \"(.+?)\";",arr[8])[0]
			ID2name_dict[gene_ID] = gene_name
		elif arr[2] == "transcript":
			gene_ID = get_right_ID(re.findall("gene_id \"(.+?)\";",arr[8])[0])
			trans_ID = get_right_ID(re.findall("transcript_id \"(.+?)\";",arr[8])[0])
			trans2gene_dict[trans_ID] = gene_ID
			
black_list_trans = []
black_list_gene = []
#### transcript to gene mapping ####
trans_prop_dict = defaultdict()
with open(isoform_prop_inf_name, 'r') as prop_inf:
	for index, line in enumerate(prop_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		gene_ID = get_right_ID(arr[1])
		trans2gene_dict[arr[0]] = gene_ID
		trans_prop_dict[arr[0]] = np.nanmedian(list(map(float, arr[2:2+tumor_count])))
		#if np.nanmedian(list(map(float, arr[2:2+tumor_count]))) < interested_isoform_prop_cutoff:
		#	black_list_trans.append(arr[0])

### filter out transcripts have high expression is tissue ###
with open(isoform_cpm_inf_name, 'r') as cpm_inf:
	for index,line in enumerate(cpm_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		tissue_array = np.array(list(map(float, arr[tumor_count+3:len(arr)])))
		if int(np.count_nonzero(tissue_array > tissue_CPM_cutoff)) > tissue_percentage_cutoff*len(tissue_array):
			black_list_trans.append(arr[0])

#### filter out genes that annotation rate is less than 50% in all tissue samples ####
with open(annotated_isoform_contri_inf_name, 'r') as gene_inf:
	for index, line in enumerate(gene_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		map_anno_rate_in_tissue = np.max(list(map(float, arr[tumor_count+1:len(arr)])))
		#if map_anno_rate_in_tissue <= anno_percent_cutoff:
		#	black_list_gene.append(arr[0])

'''
#### filter out genes that have more than 10% expression from other isoforms ####
with open(other_isoform_contri_inf_name, 'r') as gene_inf_2:
	for index, line in enumerate(gene_inf_2):
		arr = line.strip().split('\t')
		if index == 0: continue
		if float(arr[3]) >= other_isoform_prop_cutoff:
		#if np.max([float(arr[2]), float(arr[3])]) >= other_isoform_prop_cutoff:
			if arr[0] not in black_list_gene:
				black_list_gene.append(arr[0])
print ("number of black list genes:", len(black_list_gene))
'''

[n_start, n_out_1, n_out_3] = [0,0,0]
outf = open("%s/5_5_Summarized_CAR_T_prioritized_targets_temp.txt" % outf_dir, "w")
with open(tumor_vs_normal_inf_name, "r") as inf:
	for index,line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0: 
			outf.write(line)
			continue
		n_start += 1
		gene_ID = get_right_ID(arr[1])
		trans_ID = get_right_ID(arr[0])
		gene_name = ID2name_dict[gene_ID] if gene_ID in ID2name_dict else gene_ID
		## Tissue prevalence assessment ##
		if trans_ID in black_list_trans:
			n_out_1 += 1
			continue
		## Tissue prevalence assessment ##
		if arr[-3] not in ['Identified by both: consistent', 'Annotated_cell_surface:TM', 'Annotated_cell_surface:Other', 'Annotated_cell_surface:membrane_protein']: 
			n_out_3 += 1
			continue
		#if gene_ID in black_list_gene: continue		
		outf.write(line)
outf.close()

print("Number of targets before filtering\t"+str(n_start)+'\n')
print("Number of targets that have been filtered out in step1: Tissue prevalence assessment\t"+str(n_out_1)+'\n')
print("Number of targets that have been filtered out in step3: Peptide surface presentation assessment\t"+str(n_out_3)+'\n')

