import os,re,sys
from collections import defaultdict
import numpy as np

target_trans = sys.argv[1]
CPM_inf_name = sys.argv[2]
tumor_sample_count = int(sys.argv[3])
out_dir = sys.argv[4]
#target_trans = "ESPRESSO:chr7:10297:45@Melanoma"
window_size = int(sys.argv[5])
step_1_outf_name = out_dir+"/7_1_peptide_matrix_WindowSize_%s_%s.txt" % (str(window_size),target_trans)

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
		trans_ID = get_right_ID(arr[0])
		gene_ID = get_right_ID(arr[2])
		exp_dict[trans_ID] = np.array(list(map(float, arr[3:len(arr)])))

###### output sample list for each group #######
outf_sample = open(out_dir+"/7_3_sample_list.txt", 'w')
outf_sample.write('Tumor\t'+';'.join(sample_list[0:tumor_sample_count])+'\n')
outf_sample.write('Normal_tissue\t'+';'.join(sample_list[tumor_sample_count:len(sample_list)])+'\n')
outf_sample.close()
################################################


peptide2exp = defaultdict(lambda: [])
pos2peptide = defaultdict()
with open(step_1_outf_name,'r') as inf_1:
	for index,line in enumerate(inf_1):
		arr = line.strip().split('\t')
		if index == 0: continue
		this_peptide = arr[0]
		for each_pos in arr[1].split(';'):
			pos2peptide[str(each_pos)] = arr[0]
		only_exp_list = np.zeros(len(sample_list))
		for each_trans in arr[2].split(';'):
			if len(exp_dict[each_trans]) == 0: continue
			only_exp_list += exp_dict[each_trans]
		peptide2exp[this_peptide] = only_exp_list 


#start_pos = int(window_size/2)   #### now we use 1 as the start_pos
start_pos = 0   #### since we always add 1 in later steps
full_matrix = np.zeros(shape=(1,len(sample_list)))
full_matrix_log2 = np.matrix(['Peptide','Index']+sample_list)
corrected_index_list = []
for i in range(1, len(pos2peptide.keys())+1):
	this_peptide = pos2peptide[str(i)]
	only_exp_list = peptide2exp[this_peptide]
	full_matrix = np.vstack((full_matrix, only_exp_list))
	only_exp_list_log = np.log2(np.array(only_exp_list)+1)
	only_log2_list = [this_peptide,str(i+start_pos)] + list(map(str,only_exp_list_log))
	full_matrix_log2 = np.vstack((full_matrix_log2, only_log2_list))
	corrected_index_list.append(str(i+start_pos))

pseudo_count = 0.1
print ('tumor_sample', tumor_sample_count)
print ('normal sample', len(sample_list) - tumor_sample_count)
full_matrix_2 = full_matrix.transpose()
tumor_mean_row = np.mean(full_matrix_2[0:tumor_sample_count,0:], axis=0)
tissue_mean_row = np.mean(full_matrix_2[tumor_sample_count:len(sample_list),0:], axis=0)
FC_mean_row = np.log2((np.array(tumor_mean_row)+pseudo_count) / (np.array(tissue_mean_row)+pseudo_count))
#FC_mean_row_z_score = (FC_mean_row - np.mean(FC_mean_row)) / np.std(FC_mean_row)
FC_mean_index = [i for i,x in enumerate(FC_mean_row) if x == max(FC_mean_row)]
FC_mean_peptide = [pos2peptide[str(i+1)] for i in FC_mean_index]

tumor_median_row = np.median(full_matrix_2[0:tumor_sample_count,0:], axis=0)
tissue_median_row = np.median(full_matrix_2[tumor_sample_count:len(sample_list),0:], axis=0)
FC_median_row = np.log2((np.array(tumor_median_row)+pseudo_count) / (np.array(tissue_median_row)+pseudo_count))
#FC_median_row_z_score = (FC_median_row - np.mean(FC_median_row)) / np.std(FC_median_row)
FC_median_index = [i for i,x in enumerate(FC_median_row) if x == max(FC_median_row)]
FC_median_peptide = [pos2peptide[str(i+1)] for i in FC_median_index]

##### Output region with the highest tumor-specificity score #######
outf_region_name = out_dir+"/7_3_highest_score_region_WindowSize_%s_%s.txt" % (str(window_size),target_trans)
outf_region = open(outf_region_name, 'w')
print('mean', np.array(FC_mean_index)+1+start_pos, FC_mean_peptide)
print('median', np.array(FC_median_index)+1+start_pos, FC_median_peptide)
outf_region.write('Mean\t'+','.join(map(str,np.array(FC_mean_index)+1+start_pos))+'\t'+','.join(map(str,FC_mean_peptide))+'\n')
outf_region.write('Median\t'+','.join(map(str,np.array(FC_median_index)+1+start_pos))+'\t'+','.join(map(str,FC_median_peptide))+'\n')
outf_region.close()
##################

full_matrix_log_2 = full_matrix_log2.transpose()
#full_matrix_log_2_rep = full_matrix_log_2
full_matrix_log_2 = np.vstack((full_matrix_log_2, ['Log2FC_mean']+list(map(str, FC_mean_row))))
full_matrix_log_2 = np.vstack((full_matrix_log_2, ['Log2FC_median']+list(map(str, FC_median_row))))

#np.print(full_matrix_2)
outf_name = out_dir+"/7_3_calculate_score_matrix_WindowSize_%s_%s.txt" % (str(window_size),target_trans)
with open(outf_name,'w') as outf:
	for line in full_matrix_log_2:
		np.savetxt(outf, line, delimiter='\t', fmt='%s')


header_index_list = []
max_pos = 0
outf_final_name = re.sub('.txt','_reshaped.txt',outf_name)
outf_final = open(outf_final_name,'w')
outf_final.write('Index\tSample\tValue\n')
###### change format ######
with open(outf_name,'r') as inf_3:
	for index,line in enumerate(inf_3):
		arr = line.strip().split('\t')
		if index == 0: continue
		elif index == 1:
			header_index_list = arr[1:len(arr)]
			max_pos = int(arr[-1])
		else:
			this_sample = arr[0]
			for i in range(1, start_pos+1):
				outf_final.write(str(i)+'\t'+this_sample+'\tNA\n')
			for i2 in range(1,len(arr)):
				value = arr[i2]
				this_index = header_index_list[i2-1]
				outf_final.write(str(this_index)+'\t'+this_sample+'\t'+str(value)+'\n')
			for i3 in range(max_pos+1, max_pos+start_pos+1):
				outf_final.write(str(i3)+'\t'+this_sample+'\tNA\n')
outf_final.close()


last_value = 0
last_index = 0
collapsed_index = ''
collapsed_index_list = []
sample_used_list = []
col_index2value_dict = defaultdict()

####### collapsed #######

with open(outf_name,'r') as inf_4:
	for index,line in enumerate(inf_4):
		arr = line.strip().split('\t')
		if index <= 1: continue
		else:
			this_sample = arr[0]
			sample_used_list.append(this_sample)
			if this_sample != 'Log2FC_median': continue
			for i in range(1, len(arr)):
				this_value = round(float(arr[i]),4)
				this_index = header_index_list[i-1]
				if this_value != last_value:
					if collapsed_index != '':
						collapsed_index = collapsed_index+'-'+str(last_index)
						collapsed_index_list.append(collapsed_index)
						col_index2value_dict[collapsed_index] = last_value
					collapsed_index = str(this_index)
				last_index = this_index
				last_value = this_value
			#### the very last pos ###
			collapsed_index = collapsed_index+'-'+str(last_index)
			collapsed_index_list.append(collapsed_index)
			col_index2value_dict[collapsed_index] = last_value

exp_collapsed_dict = defaultdict(lambda: defaultdict(lambda: ''))
with open(outf_name,'r') as inf_5:
	for index,line in enumerate(inf_5):
		arr = line.strip().split('\t')
		if index <= 1: continue
		else:
			this_sample = arr[0]
			for each_col_index in collapsed_index_list:
				this_value = float(arr[int(each_col_index.split('-')[0])-start_pos])
				exp_collapsed_dict[this_sample][each_col_index] = this_value

######### topological annotation #####
trans2extra_pos_dict = defaultdict()
trans2tm_pos_dict = defaultdict()
with open(out_dir+'/7_2_extracellular_matrix.txt','r') as topo_inf:
	for index, line in enumerate(topo_inf):
		if index == 0: continue
		arr = line.strip().split('\t')
		trans2extra_pos_dict[arr[0]] = arr[2]
		trans2tm_pos_dict[arr[0]] = arr[1]

def decide_topo(each_region, extra_region, tm_region):
	if re.findall('-', each_region):
		[each_left, each_right] = list(map(int, each_region.split('-')))
		topo = 'Cytoplasmic'
		for each_extra_region in extra_region.split(','):
			[each_extra_left, each_extra_right] = list(map(int, each_extra_region.split('-')))
			if each_extra_left <= each_left <= each_right <= each_extra_right:
				topo = 'Extracellular'
			elif (each_extra_left <= each_left <= each_extra_right) or (each_extra_left <= each_right <= each_extra_right):
				topo = 'TM'
			else:
				continue
		if topo == 'Cytoplasmic':
			for each_tm_region in tm_region.split(','):
				if each_tm_region == 'none': continue
				[each_tm_left, each_tm_right] = list(map(int, each_tm_region.split('-')))
				if (each_tm_left <= each_left <= each_tm_right) or (each_tm_left <= each_right <= each_tm_right):
					topo = 'TM'
		return topo
	else:
		each_region = int(each_region)
		topo = 'Cytoplasmic'
		for each_extra_region in extra_region.split(','):
			[each_extra_left, each_extra_right] = list(map(int, each_extra_region.split('-')))
			if each_extra_left <= each_region <= each_extra_right:
				topo = 'Extracellular'
				break
		if topo == 'Cytoplasmic':
			for each_tm_region in tm_region.split(','):
				if each_tm_region == 'none': continue
				[each_tm_left, each_tm_right] = list(map(int, each_tm_region.split('-')))
				if each_tm_left <= each_region <= each_tm_right:
					topo = 'TM'
					break
		return topo
				
'''
######### write output #############
outf_final_name_2 = re.sub('.txt','_reshaped_collapsed_z.txt',outf_name)
outf_final_2 = open(outf_final_name_2,'w')
outf_final_2.write('Index\tSample\tValue\n')

outf_final_name_2_5 = re.sub('.txt','_collapsed_z.txt', outf_name)
outf_final_2_5 = open(outf_final_name_2_5, 'w')
outf_final_2_5.write('Sample')
for each_col_index in collapsed_index_list:
	outf_final_2_5.write('\t'+each_col_index)
outf_final_2_5.write('\n')

outf_final_2_5.write('Location')
for each_col_index in collapsed_index_list:
	topo_status = 'NA'
	if target_trans in trans2extra_pos_dict:
		topo_status = decide_topo(each_col_index, trans2extra_pos_dict[target_trans], trans2tm_pos_dict[target_trans])
	outf_final_2_5.write('\t'+topo_status)
outf_final_2_5.write('\n')


for each_sample in sample_used_list:
	outf_final_2_5.write(each_sample)
	for each_col_index in collapsed_index_list:
		row_value_list = np.array(list(map(float,exp_collapsed_dict[each_sample].values())))
		row_value_list[row_value_list > 1e308] = 0
		log2_value = float(exp_collapsed_dict[each_sample][each_col_index])
		#z_score = (float(exp_collapsed_dict[each_sample][each_col_index]) - np.nanmean(row_value_list)) / np.nanstd(row_value_list)
		outf_final_2.write(each_col_index+'\t'+each_sample+'\t'+str(log2_value)+'\n')
		outf_final_2_5.write('\t'+str(log2_value))
	outf_final_2_5.write('\n')
outf_final_2.close()
outf_final_2_5.close()
'''

##################### with topology, CPM value #############################
outf_final_name_1_5 = re.sub('.txt','_with_topology.txt', outf_name)
outf_final_1_5 = open(outf_final_name_1_5, 'w')
outf_final_1_5.write('Sample')
for each_col_index in corrected_index_list:
    outf_final_1_5.write('\t'+each_col_index)
outf_final_1_5.write('\n')

outf_final_1_5.write('Location')
for each_col_index in corrected_index_list:
	topo_status = 'NA'
	if target_trans in trans2extra_pos_dict:
		topo_status = decide_topo(each_col_index, trans2extra_pos_dict[target_trans], trans2tm_pos_dict[target_trans])
	outf_final_1_5.write('\t'+topo_status)
outf_final_1_5.write('\n')

with open(outf_name,'r') as inf_1_5:
	for line in inf_1_5:
		arr = line.strip().split('\t')
		if arr[0] in ['Peptide','Index']: continue
		outf_final_1_5.write(line)
outf_final_1_5.close()

##################### with topology, z-score #############################
outf_final_name_1_6 = re.sub('.txt','_with_topology_z_score.txt', outf_name)
outf_final_1_6 = open(outf_final_name_1_6, 'w')
outf_final_1_6.write('Sample')
for each_col_index in corrected_index_list:
	outf_final_1_6.write('\t'+each_col_index)
outf_final_1_6.write('\n')

outf_final_1_6.write('Location')
for each_col_index in corrected_index_list:
	topo_status = 'NA'
	if target_trans in trans2extra_pos_dict:
		topo_status = decide_topo(each_col_index, trans2extra_pos_dict[target_trans], trans2tm_pos_dict[target_trans])
	outf_final_1_6.write('\t'+topo_status)
outf_final_1_6.write('\n')

for each_row_index in range(2, full_matrix_log_2.shape[0]):
	each_sample = full_matrix_log_2[each_row_index,0]
	outf_final_1_6.write(each_sample)
	for each_col_index in range(1, full_matrix_log_2.shape[1]):
		col_value_list = np.array(list(map(float,full_matrix_log_2[2:int(full_matrix_log_2.shape[0])-2, each_col_index])))
		#col_value_list[col_value_list > 1e308] = 0
		log2_value = float(full_matrix_log_2[each_row_index,each_col_index])
		z_score = (log2_value - np.nanmean(col_value_list)) / np.nanstd(col_value_list)
		if each_sample in ['Log2FC_mean','Log2FC_median']:
			outf_final_1_6.write('\t'+str(log2_value))
		else:
			outf_final_1_6.write('\t'+str(z_score))
	outf_final_1_6.write('\n')
outf_final_1_6.close()

