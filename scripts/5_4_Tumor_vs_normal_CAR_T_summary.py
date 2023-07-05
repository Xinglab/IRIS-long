import os,re,sys
from collections import defaultdict
import configargparse
import numpy as np

def parse_args():
	parser = configargparse.ArgParser(description='Tumor_vs_normal_summary')
	parser.add('-ip', '--protein_inf', dest='protein_inf', type=str, help='protein_inf')
	parser.add('-ic', '--cell_surface_inf', dest='cell_surface_inf', type=str, help='cell_surface_inf')
	parser.add('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)
	parser.add('-ok', '--outf_key', dest='outf_key', type=str, help='output key word', default='summary')
	args = parser.parse_args()
	return parser, args

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return raw_ID.split('.')[0]+'-PAR-Y'
	else:
		return raw_ID.split('.')[0]

############ main script ##############
parser, args = parse_args()

outf_dir = args.outf_dir
protein_inf = args.protein_inf
cell_surface_inf = args.cell_surface_inf
key_word = args.outf_key
file_dir = os.path.dirname(os.path.realpath(__file__))

union_gene_list = []
iso_prevalence_list = []
trans2gene_dict = defaultdict()
significant_dict = defaultdict()
CPM_info_dict = defaultdict()
with open('%s/../3_2_Tumor_vs_normal_prevalence.txt' % outf_dir, 'r') as inf_1:
	for index,line in enumerate(inf_1):
		if index == 0: continue
		arr = line.strip().split('\t')
		iso_prevalence_list.append(arr[0])
		gene_ID = get_right_ID(arr[2])
		trans2gene_dict[arr[0]] = gene_ID
		if gene_ID not in union_gene_list:
			union_gene_list.append(gene_ID)
		CPM_info_dict[arr[0]] = '\t'.join([arr[0]]+[gene_ID]+arr[-6:len(arr)]+[arr[-7]]+[str(-np.log10(float(arr[-7])))])
		significant_dict[arr[0]] = 'Prevalence_test'

iso_DE_list = []
with open('%s/../3_1_Tumor_vs_normal_DE_test.txt' % outf_dir, 'r') as inf_2:
	for index,line in enumerate(inf_2):
		arr = line.strip().split('\t')
		if index == 0: 
			CPM_info_dict[arr[0]] = '\t'.join(arr)
			continue
		iso_DE_list.append(arr[0])
		gene_ID = get_right_ID(arr[1])
		arr[1] = gene_ID
		trans2gene_dict[arr[0]] = gene_ID
		if gene_ID not in union_gene_list:
			union_gene_list.append(gene_ID)
		CPM_info_dict[arr[0]] = '\t'.join(arr)
		if arr[0] in significant_dict:
			significant_dict[arr[0]] = 'Both_test'
		else:
			significant_dict[arr[0]] = 'DE_test'

########################
Union_list = list(set(iso_prevalence_list + iso_DE_list))
Intersection_list = set(iso_prevalence_list).intersection(iso_DE_list)
only_iso_prevalence = set(iso_prevalence_list).difference(iso_DE_list)
only_iso_DE = set(iso_DE_list).difference(iso_prevalence_list)
print (len(iso_prevalence_list), len(iso_DE_list))
print ('Union_list', len(Union_list), 'genes', len(union_gene_list))
print ('Intersection_list', len(Intersection_list))
print ('Only_pass_prevalence', len(only_iso_prevalence))
print ('Only_pass_DE', len(only_iso_DE))

union_out = open("%s/../3_3_Tumor_vs_normal_DE_and_prevalence_test.txt" % outf_dir, 'w')
union_out.write("Union_list\t"+str(len(Union_list))+'\t'+','.join(Union_list)+'\n')
union_out.write("Intersection_list\t"+str(len(Intersection_list))+'\t'+','.join(Intersection_list)+'\n')
union_out.write("Only_pass_prevalence\t"+str(len(only_iso_prevalence))+'\t'+','.join(only_iso_prevalence)+'\n')
union_out.write("Only_pass_DE\t"+str(len(only_iso_DE))+'\t'+','.join(only_iso_DE)+'\n')
union_out.close()
############################

Union_gene_list = []
for each_iso in Union_list:
	gene_ID = trans2gene_dict[each_iso]
	if gene_ID not in Union_gene_list:
		Union_gene_list.append(gene_ID)
print ('Union_gene_list', len(Union_gene_list))

pro_list = []
pro_gene_list = []
with open(protein_inf, 'r') as pro_inf:
	for line in pro_inf:
		if line.startswith('>'):
			trans_ID = line.split(' ')[0].split('|')[-1].split('_')[0]
			gene_ID = line.split(' ')[1]
			if trans_ID in Union_list:
				pro_list.append(trans_ID)
				if gene_ID not in pro_gene_list:
					pro_gene_list.append(gene_ID)
print ('pro_trans_list', len(pro_list))
print ('pro_gene_list', len(pro_gene_list))



outf = open('%s/5_4_Tumor_vs_normal_CAR_T_%s.txt' % (outf_dir, key_word), 'w')
outf.write(CPM_info_dict['Transcript_ID']+'\tSignificant_in\tGene_name\tdb\tUniProt_ID\tTM\tTM_seq\tExtra_pos\tExtra_topo_pos\tDescription\tType\tPredicted_topology\tConfidence\n')
surface_pro_list = []
surface_pro_gene_list = []
with open(cell_surface_inf, 'r') as surface_inf:
	for line in surface_inf:
		arr = line.strip().split('\t')
		#if arr[-1] == 'High_confidence' and arr[-3] != 'Annotated_cell_surface:Other':
			#if arr[-3] not in ['Annotated_cell_surface:TM','Identified by both: consistent'] : continue
		trans_ID = arr[0]
		gene_ID = arr[1]
		if trans_ID in pro_list:
			surface_pro_list.append(trans_ID)
			outf.write(CPM_info_dict[trans_ID]+'\t'+significant_dict[trans_ID]+'\t'+'\t'.join(arr[2:len(arr)])+'\n')
			if gene_ID not in surface_pro_gene_list:
				surface_pro_gene_list.append(gene_ID)
outf.close()
print ('surface_pro_list', len(surface_pro_list))
print ('surface_pro_gene_list', len(surface_pro_gene_list))


