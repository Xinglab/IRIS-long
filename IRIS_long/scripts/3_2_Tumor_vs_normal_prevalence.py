import os,re,sys
from collections import defaultdict
import numpy as np
import scipy
from scipy import stats
import configargparse

def parse_args():
	parser = configargparse.ArgParser(description='Generate bar and structure figures for examples')
	parser.add('-i', '--CPM_inf', dest='CPM_inf', type=str, help='CPM_inf', required=True)
	parser.add('-n', '--Tumor_num', dest='Tumor_num', type=int, help='Tumor_num', required=True)
	#parser.add('-cp', '--cutoff_p', dest='cutoff_p', type=float, default='1e-6', help='cutoff_p')
	parser.add('-ct', '--cutoff_tumor_cpm', dest='cutoff_tumor_cpm', type=float, help='cutoff_tumor_cpm')
	parser.add('-cn', '--cutoff_tissue_cpm', dest='cutoff_tissue_cpm', type=float, help='cutoff_tissue_cpm')
	parser.add('-ctp', '--cutoff_tumor_percentage', dest='cutoff_tumor_percentage', type=float, help='cutoff_tumor_percentage')
	parser.add('-cnp', '--cutoff_tissue_percentage', dest='cutoff_tissue_percentage', type=float, help='cutoff_tissue_percentage')
	parser.add('-et', '--exclude_tissue', dest='exclude_tissue', type=str, help='Excluded tissues')
	parser.add('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

CPM_inf = args.CPM_inf
Tumor_num = int(args.Tumor_num)
outf_dir = args.outf_dir
#cutoff_p = float(args.cutoff_p)
cutoff_tumor_cpm = float(args.cutoff_tumor_cpm)
cutoff_tissue_cpm = float(args.cutoff_tissue_cpm)
cutoff_tumor_percentage = float(args.cutoff_tumor_percentage)
cutoff_tissue_percentage = float(args.cutoff_tissue_percentage)
exclude_tissue_list = args.exclude_tissue.split(",")
file_dir = os.path.dirname(os.path.realpath(__file__))

####
DE_transcript_dict = defaultdict()
with open('%s/3_1_Tumor_vs_normal_DE_test.txt' % outf_dir, 'r') as inf_extra:
	for index, line in enumerate(inf_extra):
		arr = line.strip().split('\t')
		if index == 0: continue
		DE_transcript_dict[arr[0]] = 1

####
excluded_index_list = []
Tumor_sample_list = []
Tissue_sample_list = []
pseudo = 0.1
outf = open('%s/3_2_Tumor_vs_normal_prevalence.txt' % outf_dir, 'w')
with open(CPM_inf, 'r') as inf:
	for index, line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0:
			excluded_index_list = [i for i in range(len(arr)) if arr[i].title() in exclude_tissue_list]
			Tumor_sample_list = arr[3:3+Tumor_num]
			Tissue_sample_list = [arr[i] for i in range(3+Tumor_num,len(arr)) if i not in excluded_index_list] 
			outf.write(line.strip()+'\tTumor\tTissue\tP_value\tMean_Tumor\tMean_Tissue\tlog2_Mean_FC\tMedian_Tumor\tMedian_Tissue\tlog2_Median_FC\n')
			continue
		if re.findall('_PAR_Y', arr[0]):
			trans_ID = arr[0].split('.')[0]+'-PAR-Y'
		else:
			trans_ID = arr[0].split('.')[0]
		if re.findall('_PAR_Y', arr[2]):
			gene_ID = arr[2].split('.')[0]+'-PAR-Y'
		else:
			gene_ID = arr[2].split('.')[0]
		## must be tumor-enriched transcript ##
		if trans_ID not in DE_transcript_dict: continue
		Tumor_list = np.array(list(map(float,arr[3:3+Tumor_num])))
		Tumor_with_iso = len(Tumor_list[Tumor_list >= cutoff_tumor_cpm])
		Tumor_without_iso = len(Tumor_sample_list) - Tumor_with_iso
		Tumor_percentage = float(Tumor_with_iso)/len(Tumor_list)
		filtered_tissue = [arr[i] for i in range(3+Tumor_num,len(arr)) if i not in excluded_index_list]
		Tissue_list = np.array(list(map(float, filtered_tissue)))
		Tissue_with_iso = len(Tissue_list[Tissue_list >= cutoff_tissue_cpm])
		Tissue_without_iso = len(Tissue_sample_list) - Tissue_with_iso
		Tissue_percentage = float(Tissue_with_iso)/len(Tissue_list)
		table = [[Tumor_with_iso, Tissue_with_iso], [Tumor_without_iso, Tissue_without_iso]]
		p_value = stats.fisher_exact(table, alternative='two-sided')[1]
		mean_fold = np.log2((np.mean(Tumor_list)+pseudo)/(np.mean(Tissue_list)+pseudo))
		median_fold = np.log2((np.median(Tumor_list)+pseudo)/(np.median(Tissue_list)+pseudo))
		
		#if p_value < cutoff_p:
		if (Tumor_percentage >= cutoff_tumor_percentage) and (Tissue_percentage <= cutoff_tissue_percentage):
			Tumor_item = str(Tumor_with_iso)+'/'+str(len(Tumor_list))+'/'+str(round(Tumor_percentage,2))
			Tissue_item = str(Tissue_with_iso)+'/'+str(len(Tissue_list))+'/'+str(round(Tissue_percentage,2))
			outf.write(line.strip()+'\t'+Tumor_item+'\t'+Tissue_item+'\t'+str(p_value)+'\t'+str(np.mean(Tumor_list))+'\t'+str(np.mean(Tissue_list))+'\t'+str(mean_fold)+'\t'+str(np.median(Tumor_list))+'\t'+str(np.median(Tissue_list))+'\t'+str(median_fold)+'\n')
outf.close()
	


