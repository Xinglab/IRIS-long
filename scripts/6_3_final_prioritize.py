import os,re,sys
from collections import defaultdict, OrderedDict
import numpy as np

tumor_count = int(sys.argv[1])
score_cutoff = float(sys.argv[2])
tissue_CPM_cutoff = float(sys.argv[3])
tissue_number_cutoff = int(sys.argv[4])
binding_affi_cutoff = float(sys.argv[5])
isoform_proportion_inf = sys.argv[6]
annotated_isoform_contri_inf = sys.argv[7]
window_size = str(sys.argv[8])
genome = sys.argv[9]
outf_dir = sys.argv[10]
anno_percent_cutoff = 50

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
if genome in ['hg38','GRCh38']:
	gtf_inf_name = "%s/references/gencode.v39.annotation.gtf" % dir_path
elif genome in ['hg19', 'GRCh37']:
	gtf_inf_name = "%s/references/gencode.v34lift37.annotation.gtf" % dir_path
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
			
#### transcript to gene mapping ####
trans_prop_dict = defaultdict()
with open(isoform_proportion_inf, 'r') as prop_inf:
	for index, line in enumerate(prop_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		gene_ID = get_right_ID(arr[1])
		trans2gene_dict[arr[0]] = gene_ID
		trans_prop_dict[arr[0]] = np.nanmedian(list(map(float, arr[2:2+tumor_count])))

#### filter out genes that annotation rate is less than 50% in all tissue samples ####
black_list_gene = []
with open(annotated_isoform_contri_inf, 'r') as gene_inf:
	for index, line in enumerate(gene_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		map_anno_rate_in_tissue = np.max(list(map(float, arr[tumor_count+1:len(arr)])))
		if map_anno_rate_in_tissue <= anno_percent_cutoff:
			black_list_gene.append(arr[0])

#### peptide to transcripts mapping ####
peptide2HLA_dict = defaultdict()
peptide2trans_dict = defaultdict()
with open("%s/6_2_peptide_matrix_WindowSize_%s_TCR.txt" % (outf_dir, window_size), "r") as inf_s1:
	for index,line in enumerate(inf_s1):
		arr = line.strip().split('\t')
		if index == 0: continue
		peptide2HLA_dict[arr[0]] = re.sub("\*","",arr[1])
		#peptide2trans_dict[arr[0]] = arr[2].split(';')
		transcripts = arr[2].split(";")
		peptide_positions = arr[3].split(";")
		genome_positions = arr[4].split(";")
		peptide2trans_dict[arr[0]] = {"trans_ID" : transcripts, "peptide_pos": peptide_positions, "genomic_pos": genome_positions}


peptide2score_mean_dict = defaultdict()
peptide2score_median_dict = defaultdict()
with open("%s/6_2_specificity_score_matrix_WindowSize_%s_TCR.txt" % (outf_dir, window_size), 'r') as inf_s2:
	for index,line in enumerate(inf_s2):
		arr = line.strip().split('\t')
		if index == 0: continue
		if float(arr[-1]) <= score_cutoff: continue
		if float(arr[-2]) <= score_cutoff: continue
		tissue_array = np.array(list(map(float, arr[tumor_count+1:len(arr)-2])))
		if int(np.count_nonzero(tissue_array > np.log2(tissue_CPM_cutoff+1))) > tissue_number_cutoff: continue
		#if tissue_median >= np.log2(tissue_CPM_cutoff+1): continue
		peptide2score_mean_dict[arr[0]] = float(arr[-2])
		peptide2score_median_dict[arr[0]] = float(arr[-1])

pep_list = []
trans_list = []
gene_list = []
gene_name_list = []
#outf = open("./6_Summarized_TCR_prioritized_targets_score_without_limit.txt", "w")
outf = open("%s/6_3_Summarized_TCR_prioritized_targets.txt" % outf_dir, "w")
with open("%s/6_1_Summarized_TCR_netMHCpan_reshaped_out.txt" % outf_dir, "r") as inf:
	for index,line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0: 
			outf.write('\t'.join(arr)+'\tlog2FC_mean\tlog2FC_median\tDerived_transcripts\tDerived_genes\tDerived_gene_names\tPeptide_position\tGenomic_position\n')
			continue
		peptide = arr[0]
		bind_aff = np.min(list(map(float, arr[4].split(';'))))
		if re.findall('SB', arr[-1]):
			pass
		else:
			continue
		if bind_aff >= binding_affi_cutoff: continue
		if peptide not in peptide2score_median_dict: continue
		this_pep_match_trans_list = []
		this_pep_match_gene_list = []
		this_pep_match_gene_name_list = []
		this_pep_pos = []
		this_genome_pos = []
		black_gene_counts = 0
		single_peptide_dict = peptide2trans_dict[peptide]
		for i in range(len(single_peptide_dict["trans_ID"])):
			each_trans = single_peptide_dict["trans_ID"][i]
			if each_trans not in trans2gene_dict: continue
			this_pep_match_trans_list.append(each_trans)
			each_gene = trans2gene_dict[each_trans]
			if each_gene in black_list_gene: black_gene_counts += 1
			each_gene_name = ID2name_dict[each_gene] if each_gene in ID2name_dict else each_gene
			this_pep_match_gene_list.append(each_gene)
			this_pep_match_gene_name_list.append(each_gene_name)
			this_pep_pos.append(single_peptide_dict["peptide_pos"][i])
			this_genome_pos.append(single_peptide_dict["genomic_pos"][i])

		#### Derived genes cannot all be the genes that novel isoforms contribute the most of expression
		if len(this_pep_match_gene_list) == black_gene_counts: continue
		pep_list.append(peptide)

		for i in range(len(single_peptide_dict["trans_ID"])):
			each_trans = single_peptide_dict["trans_ID"][i]
			if each_trans not in trans2gene_dict: continue
			each_gene = trans2gene_dict[each_trans]
			each_gene_name = ID2name_dict[each_gene] if each_gene in ID2name_dict else each_gene
			if each_trans not in trans_list:
				trans_list.append(each_trans)
			if each_gene not in gene_list:
				gene_list.append(each_gene)
			if each_gene_name not in gene_name_list:
				gene_name_list.append(each_gene_name)

		#this_pep_match_gene_list = list(OrderedDict.fromkeys(this_pep_match_gene_list))
		#this_pep_match_gene_name_list = list(OrderedDict.fromkeys(this_pep_match_gene_name_list))
		outf.write('\t'.join(arr)+'\t'+str(peptide2score_mean_dict[peptide])+'\t'+str(peptide2score_median_dict[peptide])+'\t'+';'.join(this_pep_match_trans_list)+'\t'+';'.join(this_pep_match_gene_list)+'\t'+';'.join(this_pep_match_gene_name_list)+'\t'+';'.join(this_pep_pos)+'\t'+';'.join(this_genome_pos)+'\n')
outf.close()

outf_stats = open("%s/6_3_Summarized_TCR_prioritized_targets_stats.txt" % outf_dir, "w")
outf_stats.write("Target peptide (%sAA)\t"%window_size + str(len(pep_list))+'\n')
outf_stats.write("Target transcripts\t"+str(len(trans_list))+'\n')
outf_stats.write("Target gene\t"+str(len(gene_list))+'\n')
outf_stats.close()

