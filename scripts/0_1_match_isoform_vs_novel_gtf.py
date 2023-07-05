import os,re,sys
from collections import defaultdict


def compare(trans_1, trans_2, novel_dict, known_dict, allow_dist):  ## trans_2 is the reference transcript
	res = ['unequal','']
	gene_ID = trans2gene_dict[trans_1]
	if len(novel_dict[gene_ID][trans_1]) == len(known_dict[gene_ID][trans_2]):
		compare_flag = 0
		dist_start = 0
		dist_end = 0
		ref_isoform = ''
		len_trans_1 = 0
		len_trans_2 = 0
		strand_1 = novel_dict[gene_ID][trans_1][0].split(":")[1]
		strand_2 = novel_dict[gene_ID][trans_2][0].split(":")[1]
		if strand_1 != strand_2: 
			return ['unequal','']
		if (strand_1 == '-') and (len(novel_dict[gene_ID][trans_1]) > 1):
			if int(novel_dict[gene_ID][trans_1][0].split(':')[2]) >  int(novel_dict[gene_ID][trans_1][1].split(':')[3]):
				novel_dict[gene_ID][trans_1] = novel_dict[gene_ID][trans_1][::-1]
		if (strand_2 == '-') and (len(novel_dict[gene_ID][trans_2]) > 1):
			if int(novel_dict[gene_ID][trans_2][0].split(':')[2]) >  int(novel_dict[gene_ID][trans_2][1].split(':')[3]):
				novel_dict[gene_ID][trans_2] = novel_dict[gene_ID][trans_2][::-1]

		for i in range(len(novel_dict[gene_ID][trans_1])):
			novel_exon = ''
			known_exon = ''
			len_trans_1 += (int(novel_dict[gene_ID][trans_1][i].split(':')[3])-int(novel_dict[gene_ID][trans_1][i].split(':')[2]))
			len_trans_2 += (int(known_dict[gene_ID][trans_2][i].split(':')[3])-int(known_dict[gene_ID][trans_2][i].split(':')[2]))
			if i == 0:
				novel_dict_list = novel_dict[gene_ID][trans_1][i].split(':')
				dist_start = int(novel_dict_list[2]) - int(known_dict[gene_ID][trans_2][i].split(':')[2])
				if int(novel_dict_list[3]) == int(known_dict[gene_ID][trans_2][i].split(':')[3]):
					compare_flag += 1
			elif i == len(novel_dict[gene_ID][trans_1])-1:
				novel_dict_list = novel_dict[gene_ID][trans_1][i].split(':')
				dist_end = int(novel_dict_list[3]) - int(known_dict[gene_ID][trans_2][i].split(':')[3])
				if int(novel_dict_list[2]) == int(known_dict[gene_ID][trans_2][i].split(':')[2]):
					compare_flag += 1
			else:
				novel_exon = novel_dict[gene_ID][trans_1][i]
				known_exon = known_dict[gene_ID][trans_2][i]
				if novel_exon == known_exon:
					compare_flag += 1
		if compare_flag == len(novel_dict[gene_ID][trans_1]):
			if (abs(dist_start) <= int(allow_dist)) and (abs(dist_end) <= int(allow_dist)):
				if len_trans_1 > len_trans_2:
					ref_isoform = trans_1
					ref_len = len_trans_1
				elif len_trans_1 < len_trans_2:
					ref_isoform = trans_2
					ref_len = len_trans_2
				else:
					if index_trans_dict[trans_1] < index_trans_dict[trans_2]:
						ref_isoform = trans_1
						ref_len = len_trans_1
					else:
						ref_isoform = trans_2
						ref_len = len_trans_2
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
	identify_gtf_name = input_list[i]
	identify_abundance_name = abun_list[i]
	group = group_list[i]
	abun_inf = open(identify_abundance_name,'r')
	for index,line in enumerate(abun_inf):
		if index == 0: continue
		arr = line.strip().split('\t')
		if arr[0].startswith('ESPRESSO'):
			trans_ID = arr[0]+'@'+group
			if arr[2] == 'NA': continue
			if re.findall(',', arr[2]): continue
			gene_ID = get_right_ID(arr[2])
			trans2gene_dict[trans_ID] = gene_ID
		else:
			trans_ID = get_right_ID(arr[0])
		abundance_collect_trans_dict[group][trans_ID] = 1
	abun_inf.close()

###### load all gtfs ##########
#anno_dict = defaultdict()
identify_dict = defaultdict(lambda: defaultdict(lambda: []))
identify_trans_list = []
identify_trans_dict = defaultdict(lambda: defaultdict(lambda: ''))
#trans_chr_list_dict = defaultdict()
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
			if trans_ID.startswith('ESPRESSO'):
				trans_ID = trans_ID + '@' + group
				if trans_ID not in trans2gene_dict: continue
				gene_ID = trans2gene_dict[trans_ID]
				#if gene_ID == 'NA': continue
				this_exon = arr[0]+':'+arr[6]+':'+arr[3]+':'+arr[4]
				identify_dict[gene_ID][trans_ID].append(this_exon)
		elif arr[2] == 'transcript':
			if len(arr) > 12: 
				print (arr)
				continue
			trans_ID = get_right_ID(re.findall('transcript_id "(.+?)";',arr[8])[0])
			if trans_ID.startswith('ESPRESSO'):
				trans_ID = trans_ID + '@' + group
				if trans_ID not in trans2gene_dict: continue
				gene_ID = trans2gene_dict[trans_ID]
				#if gene_ID == 'NA': continue
				this_trans = arr[0]+':'+arr[6]+':'+arr[3]+':'+arr[4]
				index += 1
				chrom_dict[trans_ID] = arr[0]
				if trans_ID not in identify_trans_list:
					identify_trans_list.append(trans_ID)
					index_trans_dict[trans_ID] = index
				identify_trans_dict[gene_ID][trans_ID] = this_trans
			gtf_collect_trans_dict[group][trans_ID] = 1
	identify_inf.close()

####
anno_dict = identify_dict
matched_trans_dict = defaultdict()
for each_iden_trans in identify_trans_list:
	gene_ID = trans2gene_dict[each_iden_trans]
	possible_trans_list = []
	possible_trans_len_list = []
	if gene_ID in anno_dict:
		for each_anno_trans in anno_dict[gene_ID]:
			if each_iden_trans == each_anno_trans: continue
			compare_res_list = compare(each_iden_trans, each_anno_trans, identify_dict, anno_dict, allow_dist)
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
