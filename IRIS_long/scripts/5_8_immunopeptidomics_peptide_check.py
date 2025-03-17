import os,re,sys
from collections import defaultdict, OrderedDict
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

CAR_T_inf = sys.argv[1] #"./5_5_Summarized_CAR_T_prioritized_targets_final.txt"
dir_path = os.path.dirname(os.path.realpath(__file__))
obj_path = os.path.dirname(os.path.realpath(CAR_T_inf))

Final_gene_list_dict = defaultdict()
with open(CAR_T_inf,"r") as inf2:
	for index,line in enumerate(inf2):
		arr = line.strip().split("\t")
		if index == 0: continue
		Final_gene_list_dict[arr[2]] = 1

outf_gene = open(f"{obj_path}/5_6_Immunopeptidomics_hit_for_CAR_T_targets.txt","w")
with open(f"{dir_path}/references/Immunopeptidomics_all_sig_psms_count.txt","r") as Immu_inf:
	for index, line in enumerate(Immu_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_gene.write(line)
		else:
			for each_gene in arr[1].split(","):
				if each_gene in Final_gene_list_dict:
					outf_gene.write(line)	
					break
outf_gene.close()

