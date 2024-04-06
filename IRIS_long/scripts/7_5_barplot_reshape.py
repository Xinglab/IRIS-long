import os,re,sys
from collections import defaultdict
import numpy as np

target_trans = sys.argv[1]
out_dir = sys.argv[2]
window_size = sys.argv[3]
start_pos = int(sys.argv[4])
end_pos = int(sys.argv[5])
tumor_count = int(sys.argv[6])

inf_name = out_dir + '/7_3_calculate_score_matrix_WindowSize_%s_%s_with_topology.txt' % (window_size,target_trans)
outf_name = out_dir + '/7_5_barplot_WindowSize_%s_%s_reshaped.txt' % (window_size,target_trans)
outf = open(outf_name, 'w')
outf.write('Sample\tAA_index\tCPM\tlog2CPM\tGroup\n')

index_list = []
count = 0
with open(inf_name, 'r') as inf:
	for index, line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0:
			index_list = arr[1:len(arr)]
			continue
		elif arr[0] in ['Location','Log2FC_mean','Log2FC_median']: continue
		count += 1
		if count <= tumor_count:
			group = 'Tumor'
		else:
			group = 'Tissue'
		sample = arr[0]
		for i in range(1,len(arr)):
			if start_pos <= int(index_list[i-1]) <= end_pos:
				CPM = 2**float(arr[i]) - 1
				outf.write(sample+'\t'+str(index_list[i-1])+'\t'+str(CPM)+'\t'+str(arr[i])+'\t'+group+'\n')
outf.close()
