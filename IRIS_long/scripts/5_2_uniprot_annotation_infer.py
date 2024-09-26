import os,re,sys
from collections import defaultdict

abundance_inf_name = sys.argv[1] # sample_list_N2_R0_abundance.esp
fasta_inf_name = sys.argv[2] # 3_PARCB_PC.fasta
used_genome = sys.argv[3]
gencode_fasta_inf_name = sys.argv[4]
outf_dir = sys.argv[5]
key_word = sys.argv[6]

dir_path = os.path.dirname(os.path.realpath(__file__))
################################
def get_right_ID(raw_ID):
	if re.findall('_PAR_Y',raw_ID):
		new_ID = raw_ID.split('.')[0]+'-PAR-Y'
	else:
		new_ID = raw_ID.split('.')[0]
	return new_ID

# 0o24-46i59-81o96-118i138-160o202-224i236-258o273-290i312 ,  1-26,77-96,165-201,256-269
def judge_equal(tmhmm_topology,xml_extracted_array):
	Topol_arr = tmhmm_topology.split('-')
	Extra_arr = []	
	for each_part in Topol_arr:
		if re.findall('o',each_part):
			extra_piece = each_part.split('o')[0]+'-'+each_part.split('o')[1]
			Extra_arr.append(extra_piece)
	if len(Extra_arr) != len(xml_extracted_array):
		return False
	for i in range(len(Extra_arr)):
		range_1 = Extra_arr[i].split('-')
		range_2 = xml_extracted_array[i].split('-')
		####
		if int(range_1[0]) == 0: ##i=0/ since prediction always start from 0
			range_1[0] = range_2[0]
		elif i == len(Extra_arr)-1:
			range_1[1] = range_2[1]
		####
		dist_left = abs(int(range_1[0]) - int(range_2[0]))
		dist_right = abs(int(range_1[1]) - int(range_2[1]))		
		if (dist_left>20) or (dist_right>20):	
			return False
	return True

def get_domain_seq(fasta_seq,domain_str):
	domain_list = domain_str.split(':')[1].split(',')
	res = []
	for each_region in domain_list:
		each_start = int(each_region.split('-')[0])
		each_end = int(each_region.split('-')[1])
		each_seq = str(fasta_seq[each_start-1:each_end])
		res.append(each_seq)
	return ','.join(res)

def get_extra_domain_seq(ENSG_ID,fasta_seq,extra_str,extra2TM_relation):
	extra_domain_list = extra_str.split(':')[1].split(',')
	i = 0
	for each_region in extra_domain_list:
		each_start = int(each_region.split('-')[0])
		each_end = int(each_region.split('-')[1])
		each_seq = str(fasta_seq[each_start-1:each_end])
		relation_key = extra2TM_relation.split(',')[i]
		if ENSG_ID not in Gene_extra_domain_seq_dict:
			Gene_extra_domain_seq_dict[ENSG_ID] = defaultdict()
		Gene_extra_domain_seq_dict[ENSG_ID][relation_key] = each_seq
		i += 1

def check_extra_equal(ENSG_ID,inferred_extra_seq):
	res_list = []
	res = ''
	inferred_extra_seq_list = inferred_extra_seq.split(',')
	annotated_seq_list = []
	for each_relation in Gene_extra_domain_seq_dict[ENSG_ID]:
		annotated_extra_seq =  Gene_extra_domain_seq_dict[ENSG_ID][each_relation]
		annotated_seq_list.append(annotated_extra_seq)
	for i in range(len(inferred_extra_seq_list)):
		if inferred_extra_seq_list[i] in annotated_seq_list:
			res = 'TMD_'+str(i+1)+':Found'
		else:
			res = 'TMD_'+str(i+1)+':Not found'
		res_list.append(res)
	#if ENSG_ID == 'ENSG00000143630':
	#	print 'Annotated_extra_seq',Gene_extra_domain_seq_dict[ENSG_ID]
	#	print 'Inferred_extra_seq',inferred_extra_seq_list
	#	print res_list
	return ','.join(res_list)

def infer_TM_pos(fasta_seq,domain_seq):
	domain_seq_list = domain_seq.split(',')
	pos_res = []
	seq_res = []
	prev_pos = 0
	for each_domain_seq in domain_seq_list:
		displayed_pos = 'none'
		seq = 'none'
		pos_start_list = [i.start() for i in re.finditer(each_domain_seq, fasta_seq)]
		for each_start in pos_start_list:
			if each_start > prev_pos:
				each_end = int(each_start)+len(each_domain_seq)
				seq = each_domain_seq
				displayed_pos = str(int(each_start)+1)+'-'+str(each_end)
				prev_pos = each_end
				break
		pos_res.append(str(displayed_pos))
		seq_res.append(str(seq))
	final_res = [','.join(pos_res),','.join(seq_res)]
	return final_res

def infer_extra_pos(inferred_TM_pos_str,extra2TM_str,fasta_seq):
	# inferred TM: none,none,103-125,149-169
	# extra2TM_str: 1-2,3-4/ 0-1,2-3,4-END
	inferred_TM_pos_list = inferred_TM_pos_str.split(',')
	extra2TM_list = extra2TM_str.split(',')
	extra_region_list = []
	extra_seq_list = []
	IN2OUT_TM_list = []
	OUT2IN_TM_list = []
	for each_extra_region in extra2TM_list:
		IN2OUT_TM_list.append(each_extra_region.split('-')[0]) #['1','3'] or ['0','2','4']
		OUT2IN_TM_list.append(each_extra_region.split('-')[1]) #['2','4'] or ['1','3','END']

	Infer_flag = 'Yes'
	Topology_str = '0'
	Current_topology_status = ''
	for i in range(len(inferred_TM_pos_list)):     # [0,1,2,3]
		each_TM_region = inferred_TM_pos_list[i]   # when i=0, indicates the first TMD
		if each_TM_region == 'none':
			continue
		else:
			left_bound = each_TM_region.split('-')[0]
			right_bound = each_TM_region.split('-')[1]
			if str(i+1) in IN2OUT_TM_list:
				if Current_topology_status == '' or Current_topology_status == 'i':
					Topology_str = Topology_str+'i'+left_bound+'-'+right_bound
					Current_topology_status = 'o'
				else:
					Infer_flag = 'No'
					break
			elif str(i+1) in OUT2IN_TM_list:
				if Current_topology_status == '' or Current_topology_status == 'o':
					Topology_str = Topology_str + 'o'+left_bound+'-'+right_bound
					Current_topology_status = 'i'
				else:
					Infer_flag = 'No'
					break
	if Current_topology_status != '':
		Topology_str = Topology_str + Current_topology_status + str(len(fasta_seq))
	else:
		Infer_flag = 'No'  # no TMD is found			

	if Infer_flag =='No':
		extra_region_list = ['none']
		extra_seq_list = ['none']
	else:
		for each_segment in Topology_str.split('-'):
			if re.findall('o',each_segment):
				left_coor = int(each_segment.split('o')[0])
				right_coor = int(each_segment.split('o')[1])
				if right_coor == len(fasta_seq):  # the last AA
					extra_region_list.append(str(left_coor+1)+'-'+str(right_coor))
					extra_seq = fasta_seq[left_coor:right_coor]
				else:
					extra_region_list.append(str(left_coor+1)+'-'+str(right_coor-1))
					extra_seq = fasta_seq[left_coor:right_coor-1]
				extra_seq_list.append(extra_seq)
	final_res =[','.join(extra_region_list),','.join(extra_seq_list)]
	return final_res

def infer_extra_pos_original(inferred_TM_pos_str,extra2TM_str,fasta_seq):
	# inferred TM: none,none,103-125,149-169
	# extra2TM_str: 1-2,3-4,5-END
	inferred_TM_pos_list = inferred_TM_pos_str.split(',')
	extra2TM_list = extra2TM_str.split(',')
	extra_region_list = []
	extra_seq_list = []
	mapped_relation_list = []
	for each_extra2TM_relation in extra2TM_list:
		first_region = each_extra2TM_relation.split('-')[0]
		second_region = each_extra2TM_relation.split('-')[1]
		left_bound = 0
		right_bound = 0
		if first_region == '0':
			if inferred_TM_pos_list[0] == 'none':
				break
			right_bound = int(inferred_TM_pos_list[int(second_region)-1].split('-')[0])-1
			left_bound = 1
		elif second_region == 'END':
			if inferred_TM_pos_list[int(first_region)-1] == 'none':
				continue
			left_bound = int(inferred_TM_pos_list[int(first_region)-1].split('-')[1])+1
			right_bound = len(fasta_seq)
		else:
			if (inferred_TM_pos_list[int(first_region)-1] == 'none') or (inferred_TM_pos_list[int(second_region)-1]=='none'):
				if (inferred_TM_pos_list[int(first_region)-1] == 'none') and (inferred_TM_pos_list[int(second_region)-1]=='none'):
					continue
				else:
					if inferred_TM_pos_list[int(second_region)-1]=='none' and each_extra2TM_relation == extra2TM_list[-1]:
						right_bound = 'TM_loss'
						left_bound = int(inferred_TM_pos_list[int(first_region)-1].split('-')[1])+1
					else:
						break
			else:
				right_bound = int(inferred_TM_pos_list[int(second_region)-1].split('-')[0])-1
				left_bound = int(inferred_TM_pos_list[int(first_region)-1].split('-')[1])+1
		extra_region = str(left_bound)+'-'+str(right_bound)
		extra_region_list.append(extra_region)
		mapped_relation_list.append(each_extra2TM_relation)
		if right_bound == 'TM_loss':
			extra_seq = fasta_seq[int(left_bound)-1:int(left_bound)+2]+'XXX' 
		else:
			extra_seq = fasta_seq[int(left_bound)-1:int(right_bound)]
		extra_seq_list.append(extra_seq)
	final_res =[','.join(extra_region_list),','.join(extra_seq_list),','.join(mapped_relation_list)]
	return final_res

def extra_TM_relation(this_extra_str,this_TM_str):
	extra_list = this_extra_str.split(':')[1].split(',')
	TM_list = this_TM_str.split(':')[1].split(',')
	res_pos_list = []
	for i in range(len(extra_list)):
		pos = ''
		extra_s = int(extra_list[i].split('-')[0])
		extra_e = int(extra_list[i].split('-')[1])
		for j in range(len(TM_list)):
			TM_s = int(TM_list[j].split('-')[0])
			TM_e = int(TM_list[j].split('-')[1])
			if extra_e < TM_s:
				pos = str(j-1+1)+'-'+str(j+1)
				res_pos_list.append(pos)
				break
			elif TM_e < extra_s:
				if j < len(TM_list)-1:
					next_TM_s = int(TM_list[j+1].split('-')[0])
					if extra_e < next_TM_s:
						pos = str(j+1)+'-'+str(j+1+1)
						res_pos_list.append(pos)
						break
				else: # j == len(TM_list)-1
					pos = str(j+1)+'-'+str('END')
					res_pos_list.append(pos)
					break
			#last_j = j+1
	return ','.join(res_pos_list)	

###### load transcript information ########
Uni2ENST = defaultdict()
Uni2ENSG = defaultdict()
ENST2Uni = defaultdict()
ENST2ENSG = defaultdict()
ENSG2Name = defaultdict()
transcript_inf = open('%s/references/Processed_UniProt/UP000005640_9606.idmapping_coverted_2.txt' % dir_path)
for line in transcript_inf:
	arr = line.strip().split('\t')
	if len(arr) < 5: continue
	canonical_ENST = arr[1]
	Uni_ID = arr[0]
	ENST = arr[4]
	ENSG = arr[2]
	ENSG2Name[ENSG] = arr[3]
	### Uni 2 ENST    #######
	Uni2ENSG[Uni_ID] = ENSG
	if ENST != '-':
		Uni2ENST[Uni_ID] = ENST
	else:
		if canonical_ENST != '-':
			Uni2ENST[Uni_ID] = canonical_ENST
	### ENST 2 Uni #######
	if re.findall(r'-\d+',Uni_ID):  #isoform
		if (canonical_ENST != '-') and (ENST != '-'):
			if re.findall(canonical_ENST,ENST): #canonical isoform
				Uni_ID = Uni_ID.split('-')[0]
			for each_ENST_ID in ENST.split(';'):
				ENST2ENSG[each_ENST_ID] = ENSG
				if each_ENST_ID not in ENST2Uni:
					ENST2Uni[each_ENST_ID] = []
				if Uni_ID not in ENST2Uni[each_ENST_ID]:
					ENST2Uni[each_ENST_ID].append(Uni_ID)
	else: ### main isoform
		if ENST != '-':
			for each_ENST_ID in ENST.split(';'):
				ENST2ENSG[each_ENST_ID] = ENSG
				if each_ENST_ID not in ENST2Uni:
					ENST2Uni[each_ENST_ID] = []
				if Uni_ID not in ENST2Uni[each_ENST_ID]:
					ENST2Uni[each_ENST_ID].append(Uni_ID)
		elif canonical_ENST != '-':
			ENST2ENSG[canonical_ENST] = ENSG
			if canonical_ENST not in ENST2Uni:
				ENST2Uni[canonical_ENST] = []
			if Uni_ID not in ENST2Uni[canonical_ENST]:
				ENST2Uni[canonical_ENST].append(Uni_ID)			
transcript_inf.close()
############# load Novel mapped gene ############
NOVEL2ENSG = defaultdict()
inf_novel = open(abundance_inf_name)
for line in inf_novel:
	arr = line.strip().split('\t')
	if arr[2] == 'NA': continue
	if not arr[0].startswith('ENST'):
		NOVEL2ENSG[arr[0]] = get_right_ID(arr[2])
inf_novel.close()


############### feature uniprot ########
Extra_dict = defaultdict()
Cyto_dict = defaultdict()
TM_dict = defaultdict()
Protein_tag_dict = defaultdict()
status_dict = defaultdict()
uniprot_sprot_inf = open('%s/references/Processed_UniProt/uniprot_sprot_human_xml_extracellular_full.txt' % dir_path)
for line in uniprot_sprot_inf:
	arr = line.strip().split('\t')
	#Uni_ID_list = [arr[0]] + arr[1].split(';')
	if arr[0] == 'Uni_ID':
		continue
	Uni_ID = arr[0]
	Protein_tag_dict[Uni_ID] = arr[7] #Cell_surface:Other; Cell_surface:TM; Membrane_protein
	if arr[3] != '':
		for each_topol in arr[3].split(';'):
			key = each_topol.split(':')[0]
			pos = each_topol.split(':')[1]
			if key == 'Extracellular':
				if Uni_ID not in Extra_dict:
					Extra_dict[Uni_ID] = []
					status_dict[Uni_ID] = 'Reviewed'
				Extra_dict[Uni_ID].append(pos)	
			elif key == 'Cytoplasmic':
				if Uni_ID not in Cyto_dict:
					Cyto_dict[Uni_ID] = []
				Cyto_dict[Uni_ID].append(pos)
	if arr[4] != '':
		if Uni_ID not in TM_dict:
			TM_dict[Uni_ID] = arr[4].split(';')
uniprot_sprot_inf.close()

############## load fasta sequence ########################
fasta_dict = defaultdict()
inf_fasta = open(fasta_inf_name,'r')
key = ''
for line in inf_fasta:
	line = line.strip()
	if line.startswith('>'):
		key = line.split(" ")[0].split("|")[-1].split("_")[0] ## transcript_ID
		#key = '_'.join(line.lstrip('>').split('_')[0:3])
	else:
		fasta_dict[key] = line
		key = ''
inf_fasta.close()

inf_annotated_cell_surface_protein_fasta_first_time = open(gencode_fasta_inf_name,'r')
store_fasta_flag = 0
ref_ENST_ID = ''
ref_fasta_dict = defaultdict()
for line in inf_annotated_cell_surface_protein_fasta_first_time:
	line = line.strip()
	if line.startswith('>'):
		store_fasta_flag = 0
		ref_ENST_ID = ''
		if used_genome in ['hg19','HG19','GRCh37']:
			this_ENST = get_right_ID(line.split('|')[1])
		elif used_genome in ['hg38','HG38','GRCh38']:
			this_ENST = get_right_ID(line.split('|')[0].lstrip('>'))
		#this_ENST = get_right_ID(line.split(' ')[0].split('|')[2].split('_')[0])
		if this_ENST in ENST2Uni:
			for each_ID in ENST2Uni[this_ENST]:
				if each_ID in Extra_dict:
					store_fasta_flag = 1
					ref_ENST_ID = this_ENST
					break
	else:
		if (store_fasta_flag == 1) and (ref_ENST_ID != ''):
			if ref_ENST_ID not in ref_fasta_dict:
				ref_fasta_dict[ref_ENST_ID] = ''
			ref_fasta_dict[ref_ENST_ID] = ref_fasta_dict[ref_ENST_ID] + line
inf_annotated_cell_surface_protein_fasta_first_time.close()

#######################################
predicted_positive = 0
predicted_positive_true = 0
reviewed_count = 0
unreviewed_count = 0
checked_equal = 0
checked_unequal = 0
Gene_TM_dict = defaultdict()
extra2TM_dict = defaultdict()
Gene_extra_domain_seq_dict = defaultdict()

########## load annotated cell-surface protein information ####### 
ENSG_anno_status_dict = defaultdict()
inf_annotated_cell_surface_protein_fasta = open(gencode_fasta_inf_name,'r')
for line in inf_annotated_cell_surface_protein_fasta:
	line = line.strip()
	if line.startswith('>'):
		if used_genome in ['hg19','HG19','GRCh37']:
			this_ENST = get_right_ID(line.split('|')[1])
			this_ENSG = get_right_ID(line.split('|')[2])
		elif used_genome in ['hg38','HG38','GRCh38']:
			this_ENST = get_right_ID(line.split('|')[0].lstrip('>'))
			this_ENSG = get_right_ID(line.split('|')[1])
		if this_ENST in ENST2Uni:
			for each_ID in ENST2Uni[this_ENST]:
				if each_ID in Extra_dict:
					this_fasta = ref_fasta_dict[this_ENST]
					this_extra_str = 'Extra:'+','.join(Extra_dict[each_ID])
					this_TM_str = 'TM:'+','.join(TM_dict[each_ID])
					TM_seq = get_domain_seq(this_fasta,this_TM_str)
					extra2TM = extra_TM_relation(this_extra_str,this_TM_str)
					extra2TM_dict[this_ENSG] = extra2TM  #0-1,2-3,4-5,6-7,8-9,10-END
					get_extra_domain_seq(this_ENSG,this_fasta,this_extra_str,extra2TM)
					if re.findall(r'\*',TM_seq):
						continue
					Gene_TM_dict[this_ENSG] = TM_seq
					ENSG_anno_status_dict[this_ENSG] = 'Isoform has been annotated with TM'
inf_annotated_cell_surface_protein_fasta.close()
#########################################


outf_name = outf_dir+"/5_2_"+key_word+"_annotation_res.txt"
outf = open(outf_name,'w')
inf = open(fasta_inf_name,'r')
outf.write("Transcript_ID\tGene_ID\tGene_symbol\tDB\tUniProt_ID\tTM_pos\tTM_seq\tExtra_pos\tExtra_topo_pos\n")
for line in inf:
	line = line.strip()
	if line.startswith('>'):
		arr = line.split(' ')
		this_ENST = arr[0].split("|")[-1].split("_")[0]
		this_ENSG = arr[1]
		this_status = "DB"
		this_key = this_ENST
		this_fasta = fasta_dict[this_key]
		if this_key in ref_fasta_dict:
			if this_fasta != ref_fasta_dict[this_key]:
				print (f"{this_ENST} of {this_ENSG} has different protein sequence between our translation method and GENCODE.")
				continue
		this_gene_name = '-'
		if this_ENSG in ENSG2Name:
			this_gene_name = ENSG2Name[this_ENSG]

		if this_ENST in ENST2Uni:
			for each_ID in ENST2Uni[this_ENST]:  # each ID is UniProt ID
				if each_ID in Extra_dict:
					this_extra_str = 'Extra:'+','.join(Extra_dict[each_ID])
					this_TM_str = 'TM:'+','.join(TM_dict[each_ID])
					#extra_seq = get_domain_seq(this_fasta,this_extra_str)
					TM_seq = get_domain_seq(this_fasta,this_TM_str)
					extra2TM = extra_TM_relation(this_extra_str,this_TM_str)
					#print each_ID,this_extra_str,this_TM_str,extra2TM
					extra2TM_dict[this_ENSG] = extra2TM  #0-1,2-3,4-5,6-7,8-9,10-END
					get_extra_domain_seq(this_ENSG,this_fasta,this_extra_str,extra2TM)
					Gene_TM_dict[this_ENSG] = TM_seq
					outf.write(this_ENST+'\t'+this_ENSG+'\t'+this_gene_name+'\t'+this_status+'\t'+each_ID+'\t'+this_TM_str+'\t'+str(TM_seq)+'\t'+this_extra_str+'\t'+extra2TM+'\n')
				elif each_ID in Protein_tag_dict:
					this_TM_str = 'No_TM_annotation'
					TM_seq = '-'
					if each_ID in TM_dict:
						this_TM_str = 'TM:'+','.join(TM_dict[each_ID])
						TM_seq = get_domain_seq(this_fasta,this_TM_str)
					outf.write(this_ENST+'\t'+this_ENSG+'\t'+this_gene_name+'\t'+this_status+'\t'+each_ID+'\t'+this_TM_str+'\t'+str(TM_seq)+'\t'+Protein_tag_dict[each_ID]+'\t-\n')
				else:
					this_TM_str = 'No_TM_annotation'
					TM_seq = '-'
					if each_ID in TM_dict:
						this_TM_str = 'TM:'+','.join(TM_dict[each_ID])
						TM_seq = get_domain_seq(this_fasta,this_TM_str)
					outf.write(this_ENST+'\t'+this_ENSG+'\t'+this_gene_name+'\t'+this_status+'\t'+each_ID+'\t'+this_TM_str+'\t'+str(TM_seq)+'\tNo_extracellular_domain\t-\n')
		else:
			outf.write(this_ENST+'\t'+this_ENSG+'\t'+this_gene_name+'\t'+this_status+'\t-\tNo_UniProt_ID\t-\tNo_UniProt_ID\t-\n')

inf.close()
outf.close()
#print extra2TM_dict['ENSG00000162337'] 
#print extra2TM_dict['ENSG00000198794']
#print extra2TM_dict['ENSG00000149577']

################# result sort ####################
ENSG_has_ENST_dict = defaultdict()
line_dict = defaultdict()

outf_sorted_name = re.sub('res.txt','res_sorted.txt',outf_name)
inf_3 = open(outf_name,'r')
outf_3 = open(outf_sorted_name,'w')
for index, line in enumerate(inf_3):
	line = line.strip()
	arr = line.split('\t')
	if index == 0:
		outf_3.write(line+'\tDescription\n') 
		continue
	status = arr[-1]
	if arr[1] not in ENSG_has_ENST_dict:
		ENSG_has_ENST_dict[arr[1]] = []		
	if arr[0] not in ENSG_has_ENST_dict[arr[1]]:
		ENSG_has_ENST_dict[arr[1]].append(arr[0])
	if arr[1] not in ENSG_anno_status_dict:
		ENSG_anno_status_dict[arr[1]] = 'No_annotation'
	if re.findall('Extra:',arr[7]):  #Uniprot annotated extra and TM
		ENSG_anno_status_dict[arr[1]] = 'Isoform has been annotated with TM'
	elif (arr[7] == 'Cell_surface:Other') or (arr[7]=='Membrane_protein'):
		pass
	else: # Proteins that don't have TM and do not have extra annotation
		if arr[1] in Gene_TM_dict: #some isoform is annotated
			this_fasta_key = arr[0]
			#print this_fasta_key,Gene_TM_dict[arr[1]]
			infer_TM = infer_TM_pos(fasta_dict[this_fasta_key],Gene_TM_dict[arr[1]])
			arr[5] = 'Inferred_TM:'+str(infer_TM[0]) #infer_TM[0] is position, [1] is seq
			arr[6] = str(infer_TM[1])
			if list(set(arr[6].split(',')))==['none']: #no inferred TM
				arr[7] = 'Inferred_extra:none'
			else:
				#print infer_TM[0],infer_TM[1],extra2TM_dict[arr[1]]
				infer_extra = infer_extra_pos(infer_TM[0],extra2TM_dict[arr[1]],fasta_dict[this_fasta_key])
				if infer_extra[0] != 'none':
					arr[7] = 'Inferred_extra:'+str(infer_extra[0])
					#print this_fasta_key,infer_extra[1],infer_extra[2]
					arr[8] = check_extra_equal(arr[1],infer_extra[1])	
	line_dict[arr[0]+'_'+arr[1]] = '\t'.join(arr[0:len(arr)]) 
inf_3.close()

for each_gene in ENSG_has_ENST_dict.keys():
	for each_transcript in ENSG_has_ENST_dict[each_gene]:
		outf_3.write(line_dict[each_transcript+'_'+each_gene]+'\t'+ENSG_anno_status_dict[each_gene]+'\n')
outf_3.close()	
