import os,re,sys
from collections import defaultdict
import numpy as np

#'L1CAM'
required_gene_ID = sys.argv[1] 
outf_dir = sys.argv[2]
required_tumor = 'All'
required_tissue = 'All'
if len(sys.argv) > 3:
	required_tumor = sys.argv[3]
	required_tumor_list = required_tumor.split(',')
if len(sys.argv) > 4:
	required_tissue = sys.argv[4]
if required_tissue in ['ClonTech','CloneTech','clontech','clonetech']:
	required_tissue_list = ["Brain","Bladder","Whole Blood","Colon","Heart","Kidney","Liver","Lung","Ovary","Pancreas","Prostate","Muscle","Small Intestine","Spleen","Stomach","Testis","Thyroid"] 

dir_path = os.path.dirname(os.path.realpath(__file__))

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return raw_ID.split('.')[0]+'-PAR-Y'
	else:
		return raw_ID.split('.')[0]

sample2group_dict = defaultdict()
with open("/mnt/isilon/xing_lab/aspera/beazhang/gene_expression/GTEx_TCGA_all/gdc_manifest.2019-09-13.txt.map2submitterID","r") as TCGA_inf:
	for index, line in enumerate(TCGA_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		sample_ID = arr[4] #TCGA-CQ-5331-01A
		tumor_type = arr[5]
		sample2group_dict[sample_ID] = tumor_type

with open("/mnt/isilon/xing_lab/aspera/beazhang/gene_expression/GTEx_TCGA_all/sample_to_tissue_2col.txt","r") as GTEx_inf:
	for index, line in enumerate(GTEx_inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		if re.findall("Cells - ", arr[1]): continue
		sample2group_dict[arr[0]] = arr[1]

sample_list = []
Index_black_list = []
outf = open(f"{outf_dir}/{required_gene_ID}_TCGA_GTEx_gene_exp.txt", "w")
outf.write("Type\tGroup\tNormalized_count\tlog2_count\n")
with open("/mnt/isilon/xing_lab/aspera/beazhang/gene_expression/GTEx_TCGA_all/normalized_counts_hugo.txt","r") as count_inf:
	for index, line in enumerate(count_inf):
		arr = line.strip().split('\t')
		if index == 0: 
			sample_list = arr[1:len(arr)]
			for i in range(1, len(arr)):
				if arr[i].startswith("TCGA"):
					sample = '-'.join(arr[i].split(".")[0:4])
				else:
					sample = arr[i]
				if sample not in sample2group_dict:
					Index_black_list.append(i)
			continue
		else:
			gene_ID = get_right_ID(arr[0])
			if gene_ID != required_gene_ID: continue
			for i in range(1, len(arr)):
				if i in Index_black_list: continue
				sample = sample_list[i-1]
				if sample.startswith("TCGA"):
					sample = '-'.join(sample.split(".")[0:4])
					db = 'TCGA'
					group = sample2group_dict[sample]
					if required_tumor != 'All':
						if group not in required_tumor_list:
							continue
				else:
					db = 'GTEx'
					group = sample2group_dict[sample]
					if required_tissue != 'All':
						if group.split(' - ')[0] not in required_tissue_list:
							continue
				outf.write(f"{db}\t{group}\t{str(arr[i])}\t{str(np.log2(float(arr[i])+1))}\n")
outf.close()

command = f"Rscript {dir_path}/Supp_2_gene_exp_box.R {required_gene_ID} {outf_dir}"
os.system(command)

