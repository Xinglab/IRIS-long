import os,re,sys
from collections import defaultdict
import numpy as np

#### load match information ###
match_inf = open(sys.argv[1],'r')
novel2anno_dict = defaultdict()
for index, line in enumerate(match_inf):
	arr = line.strip().split('\t')
	if index == 0: continue
	if arr[0].startswith('ESPRESSO'): 
		novel2anno_dict[arr[0]] = arr[1]
match_inf.close()

#### load overlapped information ###
overlap_inf = open(sys.argv[2],'r')
all_trans_dict = defaultdict()
for line in overlap_inf:
	all_trans_dict[line.strip()] = 1
overlap_inf.close()

##### load input list file #####
input_list_inf = open(sys.argv[3],'r')
input_list = []
abun_list = []
group_list = []
for line in input_list_inf:
	input_list.append(line.strip().split('\t')[1])
	group_list.append(line.strip().split('\t')[2])
input_list_inf.close()

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return raw_ID.split('.')[0] + '-PAR-Y'
	else:
		return raw_ID.split('.')[0] 

#### load exp information ###
exp_dict = defaultdict(lambda: defaultdict())
final_sample_list = []
information_dict = defaultdict()
group2sample_num_dict = defaultdict()

for file_i in range(0,len(input_list)):
	header_list = []
	exp_inf = open(input_list[file_i],'r')
	this_group = group_list[file_i]
	for index, line in enumerate(exp_inf):
		arr = line.strip().split('\t')
		if index == 0: 
			header_list = arr[3:len(arr)]
			final_sample_list += header_list
			group2sample_num_dict[this_group] = len(arr)-3
			continue
		if arr[2] == 'NA': continue
		if re.findall(',', arr[2]): continue
		trans_ID = get_right_ID(arr[0])
		if trans_ID.startswith('ESPRESSO'):
			new_trans_ID = trans_ID+'@'+this_group
			if new_trans_ID in novel2anno_dict:
				real_trans_ID = novel2anno_dict[new_trans_ID]
			else:
				continue
		else:   ## ENST ID
			real_trans_ID = trans_ID
		information_dict[real_trans_ID] = arr[1]+'\t'+arr[2]

		if this_group not in exp_dict[real_trans_ID]:
			exp_dict[real_trans_ID][this_group] = np.array(list(map(float, arr[3:len(arr)])))
		else:
			exp_dict[real_trans_ID][this_group] += np.array(list(map(float, arr[3:len(arr)])))
	exp_inf.close()


#### output ######
outf = open(sys.argv[4], 'w')
outf.write('transcript_ID\ttranscript_name\tgene_ID\t'+'\t'.join(final_sample_list)+'\n')

for each_trans in exp_dict:
	if each_trans not in all_trans_dict: continue
	outf.write(each_trans+'\t'+information_dict[each_trans])
	for each_group in group_list:
		if each_group not in exp_dict[each_trans]:	
			exp_dict[each_trans][each_group] = [0.0] * group2sample_num_dict[each_group]
		outf.write('\t'+'\t'.join(list(map(str, exp_dict[each_trans][each_group]))))
	outf.write('\n')
outf.close()

