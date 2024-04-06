import os,re,sys
from collections import defaultdict

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
group_list = []
for line in input_list_inf:
	input_list.append(line.strip().split('\t')[0])
	group_list.append(line.strip().split('\t')[2])
input_list_inf.close()

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		right_ID = raw_ID.split('.')[0]+'-PAR-Y'
	else:
		right_ID = raw_ID.split('.')[0]
	return right_ID


#### load gtf information ###
gtf_dict = defaultdict()
annotated_flag = defaultdict()
chrom_trans_start_dict = defaultdict()
chrom_list = []

for file_i in range(0,len(input_list)):
	gtf_inf = open(input_list[file_i],'r')
	this_group = group_list[file_i]
	this_trans_ID = ''
	for line in gtf_inf:
		line = line.strip()
		if line.startswith('#'): continue
		arr = line.split('\t')
		chrom = arr[0]
		old_trans_ID = get_right_ID(re.findall('transcript_id "(.+?)"',arr[8])[0])
		if old_trans_ID.startswith('ESPRESSO'):
			trans_ID = old_trans_ID+'@'+this_group
			if trans_ID not in novel2anno_dict:  continue   # this trans has no gene information
			if trans_ID != novel2anno_dict[trans_ID]: continue  # this trans could be represented by other trans
		else:
			trans_ID = old_trans_ID
		if (arr[2] == 'transcript') and (trans_ID not in annotated_flag):
			line = re.sub(old_trans_ID, trans_ID, line)
			gtf_dict[trans_ID] = [line]
			annotated_flag[trans_ID] = 1
			if chrom not in chrom_trans_start_dict:
				chrom_trans_start_dict[chrom] = defaultdict()
				chrom_list.append(chrom)
			chrom_trans_start_dict[chrom][trans_ID] = float(arr[3])
			this_trans_ID = trans_ID
		elif (arr[2] == 'transcript') and (trans_ID in annotated_flag):
			this_trans_ID = ''
		elif (arr[2] == 'exon') and (trans_ID == this_trans_ID):
			line = re.sub(old_trans_ID, trans_ID, line)
			gtf_dict[trans_ID].append(line)
	gtf_inf.close()


#### output ######
outf = open(sys.argv[4], 'w')
chrom_list = list(set(chrom_list))
for each_chrom in chrom_list:
	sorted_trans_list = sorted(chrom_trans_start_dict[each_chrom].items(), key=lambda item:item[1])
	for each_trans_item in sorted_trans_list:
		each_trans_ID = each_trans_item[0]
		if each_trans_ID not in all_trans_dict: continue
		outf.write('\n'.join(gtf_dict[each_trans_ID])+'\n')
outf.close()

