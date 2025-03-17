import os,re,sys
from collections import defaultdict, OrderedDict
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

TCR_inf = sys.argv[1] #"./6_3_Summarized_TCR_prioritized_targets.txt"
tumor_associated_trans_inf = sys.argv[2] #"../3_1_Tumor_vs_normal_DE_test_category.txt"
dir_path = os.path.dirname(os.path.realpath(__file__))
obj_path = os.path.dirname(os.path.realpath(TCR_inf))

category_dict = defaultdict()
with open(tumor_associated_trans_inf,"r") as inf1:
	for index,line in enumerate(inf1):
		arr = line.strip().split("\t")
		if index == 0: continue
		category_dict[arr[0]] = arr[-1]

Final_trans_list_dict = defaultdict()
Final_gene_list_dict = defaultdict()
Final_junction_list_dict = defaultdict()
SJ2info_dict = defaultdict(lambda: [])
with open(TCR_inf,"r") as inf2:
	for index,line in enumerate(inf2):
		arr = line.strip().split("\t")
		if index == 0: continue
		peptide = arr[0]
		transcript_list = arr[10].split(";")
		gene_list = arr[11].split(";")
		peptide_loc_list = arr[14].split(";")
		for i in range(len(transcript_list)):
			if transcript_list[i] not in category_dict: continue
			Final_trans_list_dict[transcript_list[i]] = 1
			Final_gene_list_dict[gene_list[i]] = 1
			if re.findall("#", peptide_loc_list[i]):
				new_SJ_str = peptide_loc_list[i].split("_")[0]+"_"+str(int(peptide_loc_list[i].split("#")[0].split("-")[-1])+1)+"_"+str(int(peptide_loc_list[i].split("#")[1].split("-")[0])-1)
				Final_junction_list_dict[new_SJ_str] = 1
				SJ2info_dict[new_SJ_str].append(peptide+"#"+transcript_list[i]+"#"+gene_list[i])

outf_gene = open(f"{obj_path}/6_4_GTEx_gene_tpm_for_TCR_targets.txt","w")
with open(f"{dir_path}/references/GTEx_gene_tpm_matrix.txt","r") as GTEx_gene_inf:
	for index, line in enumerate(GTEx_gene_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_gene.write(line)
		else:
			if arr[0] in Final_gene_list_dict:
				outf_gene.write(line)
outf_gene.close()

outf_trans = open(f"{obj_path}/6_4_GTEx_transcript_tpm_for_TCR_targets.txt","w")
with open(f"{dir_path}/references/GTEx_transcript_tpm_matrix.txt","r") as GTEx_trans_inf:
	for index,line in enumerate(GTEx_trans_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_trans.write(line)
		else:
			if arr[0] in Final_trans_list_dict:
				outf_trans.write(line)
outf_trans.close()

outf_SJ = open(f"{obj_path}/6_4_GTEx_junction_tpm_for_TCR_targets.txt","w")
with open(f"{dir_path}/references/GTEx_Junction_cpm_matrix.txt","r") as GTEx_SJ_inf:
	for index,line in enumerate(GTEx_SJ_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_SJ.write(arr[0]+"\tJunction_related_info\t"+"\t".join(arr[1:len(arr)])+"\n")
		else:
			if arr[0] in Final_junction_list_dict:
				SJ_info = ";".join(SJ2info_dict[arr[0]])
				outf_SJ.write(arr[0]+"\t"+SJ_info+"\t"+"\t".join(arr[1:len(arr)])+"\n")
outf_SJ.close()
