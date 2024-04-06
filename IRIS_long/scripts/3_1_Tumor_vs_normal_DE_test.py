import os,re,sys
from collections import defaultdict
import numpy as np
import scipy
from scipy import stats
import configargparse
from statsmodels.stats.multitest import multipletests

def parse_args():
	parser = configargparse.ArgParser(description='Generate bar and structure figures for examples')
	parser.add('-i', '--CPM_inf', dest='CPM_inf', type=str, help='CPM_inf', required=True)
	parser.add('-n', '--Tumor_num', dest='Tumor_num', type=int, help='Tumor_num', required=True)
	parser.add('-cp', '--cutoff_p', dest='cutoff_p', type=float, default='0.05', help='cutoff_p')
	parser.add('-ct', '--cutoff_tumor_cpm', dest='cutoff_tumor_cpm', type=float, default='3.0', help='cutoff_tumor_cpm')
	parser.add('-cf', '--cutoff_fc', dest='cutoff_fc', type=float, default='2.0', help='cutoff_fc')
	parser.add('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)
	args = parser.parse_args()
	return parser, args

############ main script ##############
parser, args = parse_args()

CPM_inf = args.CPM_inf
Tumor_num = int(args.Tumor_num)
outf_dir = args.outf_dir
cutoff_p = float(args.cutoff_p)
cutoff_tumor_cpm = float(args.cutoff_tumor_cpm)
cutoff_fc = float(args.cutoff_fc)
file_dir = os.path.dirname(os.path.realpath(__file__))


p_value_dict = defaultdict()
outf = open('%s/3_1_Tumor_vs_normal_DE_test.txt' % outf_dir, 'w')
outf.write('Transcript_ID\tGene_ID\tMean_Tumor\tMean_Tissue\tlog2_Mean_FC\tMedian_Tumor\tMedian_Tissue\tlog2_Median_FC\tFDR_p_value\tlog10_FDR_p_value\n')
outf_all = open('%s/3_1_Tumor_vs_normal_DE_test_all.txt' % outf_dir, 'w')
outf_all.write('Transcript_ID\tGene_ID\tMean_Tumor\tMean_Tissue\tlog2_Mean_FC\tMedian_Tumor\tMedian_Tissue\tlog2_Median_FC\tFDR_p_value\tlog10_FDR_p_value\tSignificance\n')
with open(CPM_inf, 'r') as inf:
	for index,line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		trans_ID = arr[0]
		tumor_list = np.array(list(map(float, arr[3:3+Tumor_num])))
		tissue_list = np.array(list(map(float, arr[3+Tumor_num:len(arr)])))
		if np.median(tumor_list) <= cutoff_tumor_cpm: continue
		wilcoxon_test = scipy.stats.ranksums(tumor_list, tissue_list, alternative='two-sided')
		p_value_dict[trans_ID] = float(wilcoxon_test[1])

####### FDR correction #######
FDR_p_dict = defaultdict()
FDR_p_list = multipletests(pvals = list(p_value_dict.values()), alpha = cutoff_p, method="fdr_bh")[1]
for index, key_trans in enumerate(list(p_value_dict.keys())):
	FDR_p_dict[key_trans] = float(FDR_p_list[index])

#sorted_p_value = sorted(p_value_dict.items(), key=lambda x:float(x[1]))
#k_list, v_list = zip(*sorted_p_value)
#rank_dict = defaultdict()
#highest_rank = 0
#for index,key in enumerate(k_list):
#	rank_dict[key] = index+1
#	this_p_value = p_value_dict[key]*len(k_list) / (index+1)
#	if this_p_value < cutoff_p:
#		highest_rank = max(index+1, highest_rank)
#print(highest_rank, len(k_list))

pseudo = 0.1
with open(CPM_inf, 'r') as inf_2:
	for index,line in enumerate(inf_2):
		arr = line.strip().split('\t')
		if index == 0: continue
		trans_ID = arr[0]
		tumor_list = np.array(list(map(float, arr[3:3+Tumor_num])))
		tissue_list = np.array(list(map(float, arr[3+Tumor_num:len(arr)])))
		if np.median(tumor_list) < cutoff_tumor_cpm: continue
		mean_fold = np.log2((np.mean(tumor_list)+pseudo)/(np.mean(tissue_list)+pseudo))
		median_fold = np.log2((np.median(tumor_list)+pseudo)/(np.median(tissue_list)+pseudo))
		#wilcoxon_test = scipy.stats.ranksums(tumor_list, tissue_list, alternative='two-sided')
		adjust_p_value = FDR_p_dict[trans_ID]
		#adjust_p_value = wilcoxon_test[1] * len(k_list) / rank_dict[arr[0]]
		log_p = -np.log10(adjust_p_value)
		tag = "Not_significant"
		if (adjust_p_value < cutoff_p):
			if (median_fold > np.log2(cutoff_fc)):
				outf.write(arr[0]+'\t'+arr[2]+'\t'+str(np.mean(tumor_list))+'\t'+str(np.mean(tissue_list))+'\t'+str(mean_fold)+'\t'+str(np.median(tumor_list))+'\t'+str(np.median(tissue_list))+'\t'+str(median_fold)+'\t'+str(adjust_p_value)+'\t'+str(log_p)+'\n')
				tag = "Significant"
		outf_all.write(arr[0]+'\t'+arr[2]+'\t'+str(np.mean(tumor_list))+'\t'+str(np.mean(tissue_list))+'\t'+str(mean_fold)+'\t'+str(np.median(tumor_list))+'\t'+str(np.median(tissue_list))+'\t'+str(median_fold)+'\t'+str(adjust_p_value)+'\t'+str(log_p)+'\t'+tag+'\n')

outf.close()
outf_all.close()
