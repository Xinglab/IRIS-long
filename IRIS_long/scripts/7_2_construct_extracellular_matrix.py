import os,re,sys
from collections import defaultdict

tm_inf = sys.argv[1]
protein_inf_name = sys.argv[2]
out_dir = sys.argv[3]

Extra_pos_dict = defaultdict()
TM_pos_dict = defaultdict()
with open(tm_inf, 'r') as pos_inf:
	for line in pos_inf:
		arr = line.strip().split('\t')
		if arr[-1] == 'Less_confidence': continue
		if arr[-3] in ['Annotated_cell_surface:TM', 'Identified by both: consistent']:
			extra_pos = arr[7].split(':')[1].split(',')
			if re.findall('0-1', arr[8]):
				extra_pos[0] = '1-' + extra_pos[0].split('-')[1]
			Extra_pos_dict[arr[0]] = extra_pos
#			print (arr[0], arr[7], extra_pos)
			tm_pos = arr[5].split(':')[1].split(',')
			TM_pos_dict[arr[0]] = tm_pos

seq_dict = defaultdict()
trans_ID = ''
with open(protein_inf_name, 'r') as pro_seq_inf:
	for line in pro_seq_inf:
		if line.startswith('>'):
			trans_ID = line.split('|')[2].split('_')[0]
		else:
			seq_dict[trans_ID] = line.strip().upper()


outf = open(out_dir+'/7_2_extracellular_matrix.txt','w')
outf.write('Transcript_ID\tTM_pos\tExtra_pos\tExtra_peptide\n')
for each_trans in Extra_pos_dict:
	this_seq = seq_dict[each_trans]
	this_extra_pos = Extra_pos_dict[each_trans]
	this_tm_pos = TM_pos_dict[each_trans]
	tm_pep_list = []
	for each_tm in this_tm_pos:
		if each_tm == 'none': continue
		tm_pep = this_seq[int(each_tm.split('-')[0])-1:int(each_tm.split('-')[1])]
		tm_pep_list.append(tm_pep)
	extra_pep_list= []
	for each_extra in this_extra_pos:
		if each_extra == 'none': continue
		extra_pep = this_seq[int(each_extra.split('-')[0])-1:int(each_extra.split('-')[1])]
		extra_pep_list.append(extra_pep)
	outf.write(each_trans+'\t'+','.join(this_tm_pos)+'\t'+','.join(this_extra_pos)+'\t'+','.join(extra_pep_list)+'\n')
outf.close()

