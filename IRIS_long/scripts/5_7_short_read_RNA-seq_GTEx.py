import os,re,sys
from collections import defaultdict, OrderedDict
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

CAR_T_inf = sys.argv[1] #"./5_5_Summarized_CAR_T_prioritized_targets_final.txt"
tumor_associated_trans_inf = sys.argv[2] #"../3_1_Tumor_vs_normal_DE_test_category.txt"
dir_path = os.path.dirname(os.path.realpath(__file__))
obj_path = os.path.dirname(os.path.realpath(CAR_T_inf))

category_dict = defaultdict()
with open(tumor_associated_trans_inf,"r") as inf1:
	for index,line in enumerate(inf1):
		arr = line.strip().split("\t")
		if index == 0: continue
		category_dict[arr[0]] = arr[-1]

Final_trans_list_dict = defaultdict()
Final_gene_list_dict = defaultdict()
with open(CAR_T_inf,"r") as inf2:
	for index,line in enumerate(inf2):
		arr = line.strip().split("\t")
		if index == 0: continue
		Final_trans_list_dict[arr[0]] = 1
		Final_gene_list_dict[arr[1]] = 1

outf_gene = open(f"{obj_path}/5_6_GTEx_gene_tpm_for_CAR_T_targets.txt","w")
with open(f"{dir_path}/references/GTEx_gene_tpm_matrix.txt","r") as GTEx_gene_inf:
	for index, line in enumerate(GTEx_gene_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_gene.write(line)
		else:
			if arr[0] in Final_gene_list_dict:
				outf_gene.write(line)
outf_gene.close()

outf_trans = open(f"{obj_path}/5_6_GTEx_transcript_tpm_for_CAR_T_targets.txt","w")
with open(f"{dir_path}/references/GTEx_transcript_tpm_matrix.txt","r") as GTEx_trans_inf:
	for index,line in enumerate(GTEx_trans_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_trans.write(line)
		else:
			if arr[0] in Final_trans_list_dict:
				outf_trans.write(line)
outf_trans.close()

