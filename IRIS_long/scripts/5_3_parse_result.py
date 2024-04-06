import os,re,sys
from collections import defaultdict

################################
# 0o24-46i59-81o96-118i138-160o202-224i236-258o273-290i312 ,  1-26,77-96,165-201,256-269
def judge_equal(tmhmm_topology,xml_extracted_array):
	Topol_arr = tmhmm_topology.split('-')
	Extra_arr = []
	protein_total_len = int(re.split("i|o",Topol_arr[-1])[-1])
	for each_part in Topol_arr:
		if re.findall('o',each_part):
			extra_piece = each_part.split('o')[0]+'-'+each_part.split('o')[1]
			Extra_arr.append(extra_piece)
	if len(Extra_arr) != len(xml_extracted_array):
		return False
	for i in range(len(Extra_arr)):
		range_1 = list(map(int, Extra_arr[i].split('-')))
		range_2 = list(map(int, xml_extracted_array[i].split('-')))
		#if re.findall('TM_loss',xml_extracted_array[i]):
		#	continue
		####
		if int(range_1[0]) == 0: ##i=0/ since prediction always start from 0
			range_1[0] = range_2[0]
		elif int(range_1[1]) == protein_total_len: ## if the last part is extra, make the end pos of tmhmm prediction same with xml
			range_1[1] = range_2[1]
		####
		overlap_len = min(range_1[1], range_2[1]) - max(range_1[0], range_2[0]) + 1
		individual_dist = max((range_1[1]-range_1[0])+1 , (range_2[1]-range_2[0])+1)
		if overlap_len < individual_dist*0.5:
			return False
	return True
		

############ load TMHMM result #########3
outf_pre = sys.argv[1]
outf_dir = sys.argv[2]
tmhmm_dict = defaultdict()
tmhmm_line_dict = defaultdict()
protein_len_dict = defaultdict()
inf_TMHMM = open(outf_dir+'/5_1_'+outf_pre+'_tmhmm_analysis.txt')
for line in inf_TMHMM:
	arr = line.strip().split('\t')
	pro_len = arr[1].split('=')[1]
	ExpAA = arr[2].split('=')[1] #should > 18
	First60 = arr[3].split('=')[1] # the aa in TM in first 60, <20
	PredHel = arr[4].split('=')[1] #number of TM helix
	if int(PredHel)>0:
		if (float(ExpAA)>18):
			ENST_ID = arr[0].split("|")[-1].split('_')[0]
			tmhmm_dict[ENST_ID] = '1'+arr[5].split('=')[1]+pro_len
			tmhmm_line_dict[ENST_ID] = line.strip()
			protein_len_dict[ENST_ID] = pro_len
inf_TMHMM.close()

#############

count_dict = defaultdict(lambda: 0)
whether_use_in_tmhmm_dict = defaultdict()
basic_information_dict = defaultdict()
outf_1 = open(outf_dir+'/5_3_'+outf_pre+'_high_confidence.txt','w')
outf_2 = open(outf_dir+'/5_3_'+outf_pre+'_less_confidence.txt','w')
outf_stat = open(outf_dir+'/5_3_'+outf_pre+'_stats.txt','w')
inf_inferrence = open(outf_dir+'/5_2_'+outf_pre+'_annotation_res_sorted.txt')
for index,line in enumerate(inf_inferrence):
	line = line.strip()
	arr = line.split('\t')
	if index == 0:
		outf_1.write(line+'\tType\tPredicted_topology\tConfidence\n')
		outf_2.write(line+'\tType\tPredicted_topology\tConfidence\n')
		continue
	ENST_ID = arr[0]
	tag = ''
	confidence = ''
	if ENST_ID in tmhmm_dict:
		whether_use_in_tmhmm_dict[ENST_ID]='yes'
		tmhmm_str = tmhmm_dict[ENST_ID]
	else:
		tmhmm_str = '-'
	## Extra categories: Cell_surface:Other; Extra; Inferred_extra; Membrane_protein; No_extracellular_domain; No_UniProt_ID
	if re.findall('Extra:',arr[7]):
		tag = 'Annotated_cell_surface:TM'
		confidence = 'High_confidence'
	elif arr[7] == 'Cell_surface:Other':
		tag = 'Annotated_cell_surface:Other'
		confidence = 'High_confidence'
	elif arr[7] == "Membrane_protein":
		tag = 'Annotated_cell_surface:membrane_protein'
		confidence = 'High_confidence'
	elif re.findall(r'Inferred_extra:\d+',arr[7]):
		if ENST_ID in tmhmm_dict:
			xml_array = arr[7].split(':')[1].split(',')
			if judge_equal(tmhmm_str,xml_array):
				tag = 'Identified by both: consistent'
				confidence = 'High_confidence'
			else:
				tag = 'Identified by both: do not consistent'
				confidence = 'Less_confidence'
		else:
			tag = 'Only identified by inference from annotation'
			confidence = 'Less_confidence'
	else:
		basic_information_dict[arr[0]] = '\t'.join(arr[0:4])
		continue
	count_dict[tag] += 1
	if confidence == "High_confidence":
		outf_1.write(line+'\t'+tag+'\t'+tmhmm_str+'\t'+confidence+'\n')
	elif confidence == 'Less_confidence':
		outf_2.write(line+'\t'+tag+'\t'+tmhmm_str+'\t'+confidence+'\n')

inf_inferrence.close()
outf_1.close()


for each_ID in tmhmm_dict:
	#if each_ID not in whether_use_in_tmhmm_dict:
	if each_ID in basic_information_dict:
		tag = 'Only identified by TMHMM'
		confidence = 'Less_confidence'
		tmhmm_str = tmhmm_dict[each_ID]
		count_dict[tag] += 1
		outf_2.write(basic_information_dict[each_ID]+'\t-\t-\t-\t-\t-\t-\t'+tag+'\t'+tmhmm_str+'\t'+confidence+'\n')
outf_2.close()


outf_stat.write("Categories\tCount\n")
for each_cate in count_dict:
	outf_stat.write(each_cate+'\t'+str(count_dict[each_cate])+'\n')
outf_stat.close()
