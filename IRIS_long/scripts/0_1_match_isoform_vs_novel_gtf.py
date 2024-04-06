import os,re,sys
from collections import defaultdict
import numpy as np

def compare(trans_1, trans_2, novel_dict, known_dict, allow_dist):  ## trans_2 is the reference transcript
	res = ['unequal','']
	gene_ID = trans2gene_dict[trans_1]
	if len(novel_dict[trans_1]) <= len(known_dict[trans_2]):
		[compare_flag, dist_start, dist_end, ref_isoform, len_trans_1, len_trans_2] = [0, 0, 0, '', 0, 0]
		strand_1 = novel_dict[trans_1][0].split(":")[1]
		strand_2 = novel_dict[trans_2][0].split(":")[1]
		if strand_1 != strand_2: return ['unequal','']
		if (strand_1 == '-') and (len(novel_dict[trans_1]) > 1):
			if int(novel_dict[trans_1][0].split(':')[2]) > int(novel_dict[trans_1][1].split(':')[3]):
				novel_dict[trans_1] = novel_dict[trans_1][::-1]
		if (strand_2 == '-') and (len(novel_dict[trans_2]) > 1):
			if int(novel_dict[trans_2][0].split(':')[2]) > int(novel_dict[trans_2][1].split(':')[3]):
				novel_dict[trans_2] = novel_dict[trans_2][::-1]

		# first ensure the exon index of reference isoform that matches the first exon of novel isoform
		first_novel_exon = novel_dict[trans_1][0].split(':')
		ref_j = 0
		ref_len = 0
		ref_isoform = trans_2
		for j in range(len(known_dict[trans_2])):
			ref_exon_detail = known_dict[trans_2][j].split(':')
			ref_len += (int(ref_exon_detail[3])-int(ref_exon_detail[2]))+1
			if int(first_novel_exon[3]) == int(ref_exon_detail[3]):
				ref_j = j

		## Thus i-th exon of novel isoform corresponds to the i+ref_j-th exon of reference isoform
		for i in range(len(novel_dict[trans_1])):
			novel_exon_detail = novel_dict[trans_1][i].split(':')
			if (i+ref_j) >= len(known_dict[trans_2]): return ['unequal','']

			ref_exon_detail = known_dict[trans_2][i+ref_j].split(':')
			len_trans_1 += (int(novel_exon_detail[3])-int(novel_exon_detail[2]))+1
			len_trans_2 += (int(ref_exon_detail[3])-int(ref_exon_detail[2]))+1
			if i == 0:
				dist_start = np.absolute(int(novel_exon_detail[2]) - int(ref_exon_detail[2]))
				if int(novel_exon_detail[3]) == int(ref_exon_detail[3]):
					compare_flag += 1
			elif i == len(novel_dict[trans_1])-1:
				dist_end = np.absolute(int(novel_exon_detail[3]) - int(ref_exon_detail[3]))
				if int(novel_exon_detail[2]) == int(ref_exon_detail[2]):
					compare_flag += 1
			else:
				if novel_dict[trans_1][i] == known_dict[trans_2][i+ref_j]:
					compare_flag += 1
		if compare_flag == len(novel_dict[trans_1]):
			#if (abs(dist_start) <= int(allow_dist)) and (abs(dist_end) <= int(allow_dist)):
			# compare transcript length only when the exon number is equal between two isoforms
			if len(novel_dict[trans_1]) == len(known_dict[trans_2]):
				if len_trans_1 > len_trans_2:
					ref_isoform = trans_1
					ref_len = len_trans_1
			res = ['equal', ref_isoform, ref_len]
	return res


def get_right_ID(raw_ID):
	if re.findall("_PAR_Y", raw_ID):
		return raw_ID.split(".")[0]+'-PAR-Y'
	else:
		return raw_ID.split(".")[0]


##### load input list file #####
input_list_inf = open(sys.argv[1],'r')
allow_dist = int(sys.argv[2])

input_list = []
abun_list = []
group_list = []
for line in input_list_inf:
	input_list.append(line.strip().split('\t')[0])
	abun_list.append(line.strip().split('\t')[1])
	group_list.append(line.strip().split('\t')[2])
input_list_inf.close()


###### load trans2gene ########
trans2gene_dict = defaultdict()
abundance_collect_trans_dict = defaultdict(lambda: defaultdict(lambda: 1))
gtf_collect_trans_dict = defaultdict(lambda: defaultdict(lambda: 1))
for i in range(0,len(input_list)):
	[identify_gtf_name, identify_abundance_name, group] = [input_list[i], abun_list[i], group_list[i]]
	abun_inf = open(identify_abundance_name,'r')
	for index,line in enumerate(abun_inf):
		if index == 0: continue
		arr = line.strip().split('\t')
		# Novel transcript will be filtered out it can mapped to multiple genes or no gene.
		if arr[2] == 'NA': continue
		if re.findall(',', arr[2]): continue
		trans_ID = get_right_ID(arr[0])
		gene_ID = get_right_ID(arr[2])
		if arr[0].startswith('ESPRESSO'):
			trans_ID = trans_ID+'@'+group
		trans2gene_dict[trans_ID] = gene_ID
		abundance_collect_trans_dict[group][trans_ID] = 1
	abun_inf.close()

###### load all gtfs ##########
### Exon information for all isoforms
anno_dict = defaultdict(lambda: defaultdict(lambda: []))
### Exon information for all novel isoforms
identify_dict = defaultdict(lambda: defaultdict(lambda: []))
index_trans_dict = defaultdict()
chrom_dict = defaultdict()
index = 0
for i in range(0,len(input_list)):
	identify_gtf_name = input_list[i]
	group = group_list[i]
	identify_inf = open(identify_gtf_name, 'r')
	for line in identify_inf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] == 'exon':
			trans_ID = get_right_ID(re.findall('transcript_id "(.+?)";',arr[8])[0])
			this_exon = ":".join([arr[0],arr[6],arr[3],arr[4]])
			if trans_ID.startswith('ESPRESSO'):
				trans_ID = trans_ID + '@' + group
				if trans_ID not in trans2gene_dict: continue
				gene_ID = trans2gene_dict[trans_ID]
				identify_dict[gene_ID][trans_ID].append(this_exon)
			else:
				gene_ID = trans2gene_dict[trans_ID]
				identify_dict[gene_ID][trans_ID].append(this_exon)
				anno_dict[gene_ID][trans_ID].append(this_exon)
		elif arr[2] == 'transcript':
			if len(arr) > 12: 
				print (arr)
				continue
			trans_ID = get_right_ID(re.findall('transcript_id "(.+?)";',arr[8])[0])
			this_trans = ":".join([arr[0],arr[6],arr[3],arr[4]])
			if trans_ID.startswith('ESPRESSO'):
				trans_ID = trans_ID + '@' + group
				if trans_ID not in trans2gene_dict: continue
				gene_ID = trans2gene_dict[trans_ID]
				index += 1
				chrom_dict[trans_ID] = arr[0]
				if trans_ID not in index_trans_dict:
					index_trans_dict[trans_ID] = index
			gtf_collect_trans_dict[group][trans_ID] = 1
	identify_inf.close()

####
matched_trans_dict = defaultdict()
for each_iden_trans in index_trans_dict:
	gene_ID = trans2gene_dict[each_iden_trans]
	possible_trans_list = []
	possible_trans_len_list = []
	if gene_ID in anno_dict:
		for each_anno_trans in anno_dict[gene_ID]:
			if each_iden_trans == each_anno_trans: continue
			compare_res_list = compare(each_iden_trans, each_anno_trans, identify_dict[gene_ID], anno_dict[gene_ID], allow_dist)
			if compare_res_list[0]=='equal':
				possible_trans_list.append(compare_res_list[1])
				possible_trans_len_list.append(int(compare_res_list[2]))
	if len(possible_trans_list) > 0:
		max_trans_index = possible_trans_len_list.index(max(possible_trans_len_list))
		final_trans = possible_trans_list[max_trans_index]
	else:
		final_trans = each_iden_trans ## return itself
	matched_trans_dict[each_iden_trans] = final_trans
print ('matched_trans_dict', len(matched_trans_dict))



##### output file ########
outf_name = sys.argv[3]
outf = open(outf_name,'w')
outf.write('Novel_transcripts\tMatched_longer_transcripts\n')
for each_trans in matched_trans_dict:
	matched_trans_ID = matched_trans_dict[each_trans]
	outf.write(each_trans+'\t'+matched_trans_ID+'\n')
outf.close()

outf2_name = sys.argv[4]
outf2 = open(outf2_name,'w')
all_list = []
for each_group in abundance_collect_trans_dict:
	seq_temp = set(abundance_collect_trans_dict[each_group].keys()).intersection(set(gtf_collect_trans_dict[each_group].keys()))
	all_list.append(seq_temp)
union_list = list(set.union(*all_list))
for each_trans in union_list:
	outf2.write(each_trans+'\n')
outf2.close()
