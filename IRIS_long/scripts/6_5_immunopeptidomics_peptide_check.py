import os,re,sys
from collections import defaultdict, OrderedDict
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

TCR_inf = sys.argv[1] #"./6_3_Summarized_TCR_prioritized_targets.txt"
dir_path = os.path.dirname(os.path.realpath(__file__))
obj_path = os.path.dirname(os.path.realpath(TCR_inf))

Final_peptide_list_dict = defaultdict()
with open(TCR_inf,"r") as inf2:
	for index,line in enumerate(inf2):
		arr = line.strip().split("\t")
		if index == 0: continue
		peptide = arr[0]
		transcript_list = arr[10].split(";")
		gene_list = arr[11].split(";")
		peptide_loc_list = arr[14].split(";")
		Final_peptide_list_dict[peptide] = 1

outf_pep = open(f"{obj_path}/6_4_Immunopeptidomics_hit_for_TCR_targets.txt","w")
with open(f"{dir_path}/references/Immunopeptidomics_all_sig_psms_count.txt","r") as Immu_inf:
	for index, line in enumerate(Immu_inf):
		arr = line.strip().split("\t")
		if index == 0:
			outf_pep.write(line)
		else:
			if arr[0] in Final_peptide_list_dict:
				outf_pep.write(line)
outf_pep.close()
