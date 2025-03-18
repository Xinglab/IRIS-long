import os,re,sys
from collections import defaultdict
import numpy as np
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests

def get_right_ID(raw_ID):
	if re.findall("_PAR_Y",raw_ID):
		return raw_ID.split(".")[0]+"-PAR-Y"
	else:
		return raw_ID.split(".")[0]

Tumor_vs_normal_DE_test_inf = sys.argv[1]
Tumor_vs_normal_prevalence_inf = sys.argv[2]
Tumor_num = int(sys.argv[3])
samples_abundance_combined_CPM_proportion_inf = sys.argv[4]
samples_abundance_combined_CPM_gene_inf = sys.argv[5]
obj_path = os.path.dirname(os.path.realpath(Tumor_vs_normal_DE_test_inf))

cutoff_p = 0.05
te_trans_dict = defaultdict(lambda: 0)
te_genes_dict = defaultdict(lambda: 0)
FDR_p_dict_trans = defaultdict()
with open(Tumor_vs_normal_DE_test_inf,"r") as te_inf:
	for index,line in enumerate(te_inf):
		arr = line.strip().split("\t")
		if index == 0: continue
		te_trans_dict[arr[0]] += 1
		te_genes_dict[get_right_ID(arr[1])] += 1
		FDR_p_dict_trans[arr[0]] = float(arr[-2])

with open(Tumor_vs_normal_prevalence_inf,"r") as ts_inf:
	for index,line in enumerate(ts_inf):
		arr = line.strip().split("\t")
		if index == 0: continue
		te_trans_dict[arr[0]] += 1
		te_genes_dict[get_right_ID(arr[2])] += 1

p_value_dict_trans = defaultdict()
proportion_dict_trans = defaultdict()
with open(samples_abundance_combined_CPM_proportion_inf,"r") as propor_inf:
	for index,line in enumerate(propor_inf):
		arr = line.strip().split("\t")
		if index == 0: continue
		trans_ID = arr[0]
		if trans_ID not in te_trans_dict: continue
		tumor_list = np.array(list(map(float, arr[2:2+Tumor_num])))
		tissue_list = np.array(list(map(float, arr[2+Tumor_num:len(arr)])))
		proportion_dict_trans[trans_ID] = np.median(tumor_list) - np.median(tissue_list)
		#wilcoxon_test = scipy.stats.ranksums(tumor_list, tissue_list, alternative='two-sided')
		#p_value_dict_trans[trans_ID] = float(wilcoxon_test[1])

######### Gene exp test ##############
p_value_dict_genes = defaultdict()
fc_dict_genes = defaultdict()
pseudo = 0.1
with open(samples_abundance_combined_CPM_gene_inf,"r") as gene_inf:
	for index,line in enumerate(gene_inf):
		arr = line.strip().split("\t")
		if index == 0: continue
		gene_ID = arr[0]
		if gene_ID not in te_genes_dict: continue
		tumor_list_gene = np.array(list(map(float, arr[1:1+Tumor_num])))
		tissue_list_gene = np.array(list(map(float, arr[1+Tumor_num:len(arr)])))
		fc_dict_genes[gene_ID] = np.log2((np.median(tumor_list_gene)+pseudo)/(np.median(tissue_list_gene)+pseudo))
		wilcoxon_test_gene = scipy.stats.ranksums(tumor_list_gene, tissue_list_gene, alternative='two-sided')
		p_value_dict_genes[gene_ID] = float(wilcoxon_test_gene[1])

####### Gene test FDR correction #######
FDR_p_dict_genes = defaultdict()
FDR_p_list_genes = multipletests(pvals = list(p_value_dict_genes.values()), alpha = cutoff_p, method="fdr_bh")[1]
for index, key_gene in enumerate(list(p_value_dict_genes.keys())):
	FDR_p_dict_genes[key_gene] = float(FDR_p_list_genes[index])

tag_count_dict = defaultdict(lambda: 0)
### Categorizing ####
outf_0 = open(f"{obj_path}/3_1_Tumor_vs_normal_DE_test_category.txt","w")
with open(Tumor_vs_normal_DE_test_inf,"r") as te_inf2:
	for index,line in enumerate(te_inf2):
		arr = line.strip().split("\t")
		if index == 0: 
			outf_0.write(line.strip()+'\tGene_FC\tGene_FDR\tProportion_diff\tCategory\n')
			continue
		trans_ID = arr[0]
		gene_ID = get_right_ID(arr[1])
		tag_list = []
		## gene expression ##
		if (FDR_p_dict_genes[gene_ID] < 0.05) and (fc_dict_genes[gene_ID] > 1):
			tag_list.append("Gene_upregulation")
		## trans proportion ##
		if proportion_dict_trans[trans_ID] > 10:
			tag_list.append("Isoform_switch")
		if len(tag_list) > 1:
			overall_tag = "Dual_classification"
		elif len(tag_list) == 1:
			overall_tag = tag_list[0]
		else:
			overall_tag = "Others"
		tag_count_dict[overall_tag] += 1
		str_list = arr + [str(fc_dict_genes[gene_ID]),str(FDR_p_dict_genes[gene_ID]),str(proportion_dict_trans[trans_ID]), overall_tag]
		outf_0.write("\t".join(str_list)+"\n")

with open(Tumor_vs_normal_prevalence_inf,"r") as ts_inf2:
	for index,line in enumerate(ts_inf2):
		arr = line.strip().split("\t")
		if index == 0:
			continue
		trans_ID = arr[0]
		gene_ID = get_right_ID(arr[2])
		if te_trans_dict[trans_ID] > 1: continue
		tag_list = []
		## gene expression ##
		if (FDR_p_dict_genes[gene_ID] < 0.05) and (fc_dict_genes[gene_ID] > 1):
			tag_list.append("Gene_upregulation")
		## trans proportion ##
		if proportion_dict_trans[trans_ID] > 10:
			tag_list.append("Isoform_switch")
		if len(tag_list) > 1:
			overall_tag = "Dual_classification"
		elif len(tag_list) == 1:
			overall_tag = tag_list[0]
		else:
			overall_tag = "Others"
		tag_count_dict[overall_tag] += 1
		str_list = [trans_ID, gene_ID] + arr[-6:] + ['-','-', str(fc_dict_genes[gene_ID]),str(FDR_p_dict_genes[gene_ID]), str(proportion_dict_trans[trans_ID]), overall_tag] 
		outf_0.write("\t".join(str_list)+"\n")
outf_0.close()

outf = open(f"{obj_path}/TE_trans_classification.txt","w")
outf.write("Category\tCount\n")
for each_tag in ["Dual_classification","Isoform_switch","Gene_upregulation","Others"]:
	outf.write(f"{each_tag}\t{tag_count_dict[each_tag]}\n")
outf.close()

