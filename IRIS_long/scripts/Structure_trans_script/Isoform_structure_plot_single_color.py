import re,os,sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as mPath
import matplotlib.patches as mpatches
import numpy as np


##### pre-set parameters ##############
fig = plt.figure(figsize=(12, 4),dpi=300)
ax = fig.add_subplot(111)

exon_height = 6
line_space = 13
font_size = 12
# parameter: plot size: 1000 * 600
####### input parameters ##################
CPM_dict = defaultdict()
exon_dict = defaultdict(lambda: [])
intron_dict = defaultdict(lambda: [])
corrected_exon_dict = defaultdict(lambda: [])
corrected_intron_dict = defaultdict(lambda: [])
corrected_CDS_exon_dict = defaultdict(lambda: [])
corrected_CDS_intron_dict = defaultdict(lambda: [])

min_pos = 10000000000
max_pos = 0
num_trans = 0
trans_ID_list = []
top_rank_isoform = []
bed_list_name = sys.argv[1]
gene_name = sys.argv[2]
gene_column = int(sys.argv[3])
CDS_inf_name = sys.argv[5]

tool = 'ESPRESSO'
strand = '+'
intron_shrinkage_fold = 10
if len(sys.argv)>4:
	tool = sys.argv[4]
if len(sys.argv)>6:
	intron_shrinkage_fold = int(sys.argv[6])

corrected_max_pos = 0
text_x = ''

def decide_strand(input_strand):
	if input_strand == '+':
		left_end = "5'"
		right_end = "3'"
	elif input_strand == '-':
		left_end = "3'"
		right_end = "5'"
	return [left_end, right_end]

def find_initial_CDS_pos(CDS_initial_pos):
	for each_exon in map_exon_2_corrected_exon_dict:
		[each_exon_original_left, each_exon_original_right] = map(float, each_exon.split(":"))
		[each_exon_corrected_left, each_exon_corrected_right] = map(float, map_exon_2_corrected_exon_dict[each_exon].split(":"))
		if each_exon_original_left <= CDS_initial_pos <= each_exon_original_right:
			CDS_corrected_pos = each_exon_corrected_left + (CDS_initial_pos - each_exon_original_left)
			return CDS_corrected_pos

strand_dict = defaultdict()
tag_dict = defaultdict()
abundance_index_dict = defaultdict()
## try to get min_pos, max_pos, corrected_max_pos and exon_dict
with open(bed_list_name,'r') as list_inf:
	for list_index,list_line in enumerate(list_inf):
		num_trans += 1
		first_arr = list_line.strip().split('\t')
		file_name = first_arr[1].split('/')[-1]
		this_strand = first_arr[2]
		file_name_list = file_name.split('.')[0].split('_')
		group = '_'.join(file_name_list[0:gene_column])
		#gene_name = file_name_list[2]
		if tool == 'ESPRESSO':
			trans_ID = '_'.join(file_name_list[gene_column+1:len(file_name_list)])
		elif tool == 'FLAMES':
			trans_ID = re.findall('Target.+ENS.+(ENS.+)',file_name.split('.')[0])[0]
		trans_ID_list.append(trans_ID)
		strand_dict[trans_ID] = this_strand
		abundance_index_dict[trans_ID] = int(first_arr[0])
		if len(first_arr)>3:
			tag_dict[trans_ID] = first_arr[3]
		with open(first_arr[1],'r') as inf:
			for index,line in enumerate(inf):
				arr = line.strip().split('\t')
				if index == 0: continue
				CPM_dict[trans_ID] = float(arr[3])
				exon_dict[trans_ID].append('#'.join(arr[0:3]))
				min_pos = min(min_pos, int(arr[1]))
				max_pos = max(max_pos, int(arr[2]))
		for i in range(len(exon_dict[trans_ID])-1):
			intron_left = int(exon_dict[trans_ID][i].split("#")[2])+1
			intron_right = int(exon_dict[trans_ID][i+1].split("#")[1])-1+1
			intron_dict[trans_ID].append([intron_left,intron_right])

## obtain the intronic region shared by all transcripts ##
shared_intron_list = intron_dict[trans_ID_list[0]]
shared_intron_set = set.union(*[set(range(shared_intron_list[i][0],shared_intron_list[i][1]+1)) for i in range(len(shared_intron_list))])
for trans_index in range(1,len(trans_ID_list)):
	cmp_intron_list = intron_dict[trans_ID_list[trans_index]]
	cmp_intron_set = set.union(*[set(range(cmp_intron_list[i][0],cmp_intron_list[i][1]+1)) for i in range(len(cmp_intron_list))])
	shared_intron_set = shared_intron_set.intersection(cmp_intron_set)

## store CDS information ##
CDS_dict = defaultdict()
with open(CDS_inf_name, 'r') as CDS_inf:
	for index,line in enumerate(CDS_inf):
		arr = line.strip().split("\t")
		if index == 0: continue
		if re.findall("PC|NMD", arr[-1]):
			trans_ID = arr[0]
			if trans_ID in abundance_index_dict:
				CDS_dict[trans_ID] = arr[5].split(";")

### sort color list ###
index2color_list_original = ['#FF5F42','#003E7F','#0068AF','#5495E1','#A1C2E8']
index2color_list_original = index2color_list_original[0:num_trans]
index2color_list = index2color_list_original[::-1]
map_exon_2_corrected_exon_dict = defaultdict()

## obtain the corrected maximal position
with open(bed_list_name,'r') as list_inf_0:
	for list_index,list_line in enumerate(list_inf_0):
		file_name = list_line.strip().split('\t')[1].split('/')[-1]
		file_name_list = file_name.split('.')[0].split('_')
		group = '_'.join(file_name_list[0:gene_column])
		#gene_name = file_name_list[2]
		if tool == 'ESPRESSO':
			trans_ID = '_'.join(file_name_list[gene_column+1:len(file_name_list)])
		elif tool == 'FLAMES':
			trans_ID = re.findall('Target.+ENS.+(ENS.+)',file_name.split('.')[0])[0]
		### adjust initial current_pos
		initial_pos = (float(exon_dict[trans_ID][0].split('#')[1]) - min_pos) + 1
		current_pos = initial_pos
		if len(exon_dict[trans_ID]) >= 2:
			for i in range(len(exon_dict[trans_ID])-1):
				exon_len = (float(exon_dict[trans_ID][i].split('#')[2]) - float(exon_dict[trans_ID][i].split('#')[1]))
				current_pos = current_pos + exon_len
				this_intron_left = int(exon_dict[trans_ID][i].split('#')[2]) + 1
				this_intron_right = int(exon_dict[trans_ID][i+1].split('#')[1])
				this_intron_set = set(range(this_intron_left, this_intron_right+1))
				shrinkage_len = len(this_intron_set.intersection(shared_intron_set))
				intron_len = (len(this_intron_set)-shrinkage_len) + float(shrinkage_len)/intron_shrinkage_fold
				current_pos = current_pos + intron_len
			last_exon_len = float(exon_dict[trans_ID][-1].split('#')[2]) - float(exon_dict[trans_ID][-1].split('#')[1])
			current_pos = current_pos+last_exon_len
		#### single exon #####
		else:
			exon_len = (float(exon_dict[trans_ID][0].split('#')[2]) - float(exon_dict[trans_ID][0].split('#')[1]))
			current_pos = current_pos + exon_len
		corrected_max_pos = max(corrected_max_pos,current_pos)

### the coordinate format in bed file is 0-based, e.g. (0,100], which represent from 1 to 100, in total 100 nts.
### Example [[0,100],[200,300]], the exons are actually [1...100],[201..300], the intron is [101..200]
with open(bed_list_name,'r') as list_inf_2:
	for list_index,list_line in enumerate(list_inf_2):
		file_name = list_line.strip().split('\t')[1].split('/')[-1]
		file_name_list = file_name.split('.')[0].split('_')
		group = '_'.join(file_name_list[0:gene_column])
		#gene_name = file_name_list[2]
		if tool == 'ESPRESSO':
			trans_ID = '_'.join(file_name_list[gene_column+1:len(file_name_list)])
		elif tool == 'FLAMES':
			trans_ID = re.findall('Target.+ENS.+(ENS.+)',file_name.split('.')[0])[0]
		### adjust initial current_pos
		initial_pos = (float(exon_dict[trans_ID][0].split('#')[1]) - min_pos) + 1
		current_pos = initial_pos
		if len(exon_dict[trans_ID]) >= 2:
			for i in range(len(exon_dict[trans_ID])-1):
				exon_len = (float(exon_dict[trans_ID][i].split('#')[2]) - float(exon_dict[trans_ID][i].split('#')[1]))
				current_exon_pos = str(current_pos)+':'+str(current_pos+exon_len-1)
				current_pos = current_pos + exon_len
				original_exon_pos = str(int(exon_dict[trans_ID][i].split('#')[1])+1)+":"+str(exon_dict[trans_ID][i].split('#')[2])
				map_exon_2_corrected_exon_dict[original_exon_pos] = current_exon_pos

				this_intron_left = int(exon_dict[trans_ID][i].split('#')[2]) + 1
				this_intron_right = int(exon_dict[trans_ID][i+1].split('#')[1])
				this_intron_set = set(range(this_intron_left, this_intron_right+1))
				shrinkage_len = len(this_intron_set.intersection(shared_intron_set))
				intron_len = (len(this_intron_set)-shrinkage_len) + float(shrinkage_len)/intron_shrinkage_fold
				current_intron_pos = str(current_pos)+':'+str(current_pos+intron_len-1)
				current_pos = current_pos + intron_len

				corrected_exon_dict[trans_ID].append(current_exon_pos)
				corrected_intron_dict[trans_ID].append(current_intron_pos)
			last_exon_len = float(exon_dict[trans_ID][-1].split('#')[2]) - float(exon_dict[trans_ID][-1].split('#')[1])
			current_exon_pos = str(current_pos)+':'+str(current_pos+last_exon_len-1)
			current_pos = current_pos+last_exon_len
			corrected_exon_dict[trans_ID].append(current_exon_pos)
			original_exon_pos = str(int(exon_dict[trans_ID][-1].split('#')[1])+1)+":"+str(exon_dict[trans_ID][-1].split('#')[2])
			map_exon_2_corrected_exon_dict[original_exon_pos] = current_exon_pos
		#### single exon #####
		else:
			exon_len = (float(exon_dict[trans_ID][0].split('#')[2]) - float(exon_dict[trans_ID][0].split('#')[1]))
			current_exon_pos = str(current_pos)+':'+str(current_pos+exon_len-1)
			current_pos = current_pos + exon_len
			corrected_exon_dict[trans_ID].append(current_exon_pos)
			original_exon_pos = str(int(exon_dict[trans_ID][0].split('#')[1])+1)+":"+str(exon_dict[trans_ID][0].split('#')[2])
			map_exon_2_corrected_exon_dict[original_exon_pos] = current_exon_pos

		#### correct CDS coordinate ####
		### CDS file use 1-based coordinate, so the coordinate of exon left should -1
		### Example [[1,100],[201,300]], the exons are actually [1...100],[201..300], the intron is [101..200]
		if trans_ID in CDS_dict:
			initial_CDS_pos = find_initial_CDS_pos(int(CDS_dict[trans_ID][0].split('#')[0]))
			current_CDS_pos = initial_CDS_pos
			if len(CDS_dict[trans_ID]) >= 2:
				for i in range(len(CDS_dict[trans_ID])-1):
					CDS_exon_len = (float(CDS_dict[trans_ID][i].split('#')[1]) - float(CDS_dict[trans_ID][i].split('#')[0]) + 1)
					current_CDS_exon_pos = str(current_CDS_pos)+':'+str(current_CDS_pos+CDS_exon_len-1)
					current_CDS_pos = current_CDS_pos + CDS_exon_len

					this_CDS_intron_left = int(CDS_dict[trans_ID][i].split('#')[1]) + 1
					this_CDS_intron_right = int(CDS_dict[trans_ID][i+1].split('#')[0]) - 1
					this_CDS_intron_set = set(range(this_CDS_intron_left, this_CDS_intron_right+1))
					CDS_shrinkage_len = len(this_CDS_intron_set.intersection(shared_intron_set))
					CDS_intron_len = (len(this_CDS_intron_set)-CDS_shrinkage_len) + float(CDS_shrinkage_len)/intron_shrinkage_fold
					current_CDS_intron_pos = str(current_CDS_pos)+':'+str(current_CDS_pos+CDS_intron_len-1)
					current_CDS_pos = current_CDS_pos + CDS_intron_len

					corrected_CDS_exon_dict[trans_ID].append(current_CDS_exon_pos)
					corrected_CDS_intron_dict[trans_ID].append(current_CDS_intron_pos)
				last_CDS_exon_len = float(CDS_dict[trans_ID][-1].split('#')[1]) - float(CDS_dict[trans_ID][-1].split('#')[0]) + 1
				current_CDS_exon_pos = str(current_CDS_pos)+':'+str(current_CDS_pos+last_CDS_exon_len-1)
				corrected_CDS_exon_dict[trans_ID].append(current_CDS_exon_pos)
				current_CDS_pos = current_CDS_pos+last_CDS_exon_len
			#### single exon #####
			else:
				CDS_exon_len = (float(CDS_dict[trans_ID][0].split('#')[1]) - float(CDS_dict[trans_ID][0].split('#')[0])) + 1
				current_CDS_exon_pos = str(current_CDS_pos)+':'+str(current_CDS_pos+CDS_exon_len-1)
				current_CDS_pos = current_CDS_pos + CDS_exon_len
				corrected_CDS_exon_dict[trans_ID].append(current_CDS_exon_pos)		

		####
		if text_x == '':
			text_x = -0.6*corrected_max_pos
		#print trans_ID, corrected_exon_dict[trans_ID], corrected_intron_dict[trans_ID], current_pos

		######## plot ######
		selected_color = index2color_list[list_index]
		#if (trans_ID.startswith('ESPRESSO') or (re.findall('_',trans_ID))):
		#	selected_color = '#E04036'
		#else:
		#	selected_color = '#2893D4'
		######## draw intron line ########	
		if len(corrected_exon_dict[trans_ID]) >= 2:
			intron_left = float(corrected_exon_dict[trans_ID][0].split(':')[1])
			intron_right = float(corrected_exon_dict[trans_ID][-1].split(':')[0])
			intron_x_region = np.linspace(intron_left, intron_right, round(max(10,0.1*(intron_right-intron_left))))
			intron_y_region = np.array([line_space*(list_index+1)]*len(intron_x_region))
			#print intron_x_region
			if strand_dict[trans_ID] == '-':
				plt.plot(-intron_x_region,intron_y_region,c="black",lw=0.5)
			else:
				plt.plot(intron_x_region,intron_y_region,c="black",lw=0.5)

		######## draw exon box ########
		for each_exon in corrected_exon_dict[trans_ID]:
			exon_left = float(each_exon.split(':')[0])
			exon_right = float(each_exon.split(':')[1])
			exon_length = exon_right - exon_left + 1
			#left_top_point_y = line_space*(list_index+1) + exon_height/2
			#exon_box = mpatches.Rectangle((exon_left,left_top_point_y),exon_length,exon_height,linewidth=2,edgecolor='r',facecolor='none')
			#ax.add_patch(exon_box)
			exon_x_region = np.linspace(exon_left, exon_right, round(max(10,0.1*exon_length)))
			exon_y_region_top = np.array([line_space*(list_index+1)+exon_height/6]*len(exon_x_region))
			exon_y_region_bottom = np.array([line_space*(list_index+1)-exon_height/6]*len(exon_x_region))
			#print exon_left,exon_right,exon_length,exon_x_region
			if strand_dict[trans_ID] == '-':
				ax.fill_between(-exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="black", linewidth=0.3, zorder=2)
			else:
				ax.fill_between(exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="black", linewidth=0.3, zorder=2)

		## only draw CDS region in the fist transcript ##
		if list_index == num_trans-1:
			if trans_ID in CDS_dict:
				for each_exon in corrected_CDS_exon_dict[trans_ID]:
					exon_left = float(each_exon.split(':')[0])
					exon_right = float(each_exon.split(':')[1])
					exon_length = exon_right - exon_left + 1
					CDS_exon_x_region = np.linspace(exon_left, exon_right, round(max(10,0.1*exon_length)))
					CDS_exon_y_region_top = np.array([line_space*(list_index+1)+exon_height/2]*len(CDS_exon_x_region))
					CDS_exon_y_region_bottom = np.array([line_space*(list_index+1)-exon_height/2]*len(CDS_exon_x_region))
					#print exon_left,exon_right,exon_length,exon_x_region
					if strand_dict[trans_ID] == '-':
						ax.fill_between(-CDS_exon_x_region, CDS_exon_y_region_top, CDS_exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="black", linewidth=0.3, zorder=3)
					else:
						ax.fill_between(CDS_exon_x_region, CDS_exon_y_region_top, CDS_exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="black", linewidth=0.3, zorder=3)

			## remove exon border ##
			for each_exon in corrected_exon_dict[trans_ID]:
				exon_left = float(each_exon.split(':')[0])
				exon_right = float(each_exon.split(':')[1])
				exon_length = exon_right - exon_left + 1
				basic_unit_to_display = corrected_max_pos / 1000
				exon_x_region = np.linspace(exon_left+basic_unit_to_display, exon_right-basic_unit_to_display, round(max(10,0.1*exon_length)))
				exon_y_region_top = np.array([line_space*(list_index+1)+exon_height/6-0.1]*len(exon_x_region))
				exon_y_region_bottom = np.array([line_space*(list_index+1)-exon_height/6+0.1]*len(exon_x_region))
				#print exon_left,exon_right,exon_length,exon_x_region
				if strand_dict[trans_ID] == '-':
					ax.fill_between(-exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="none", zorder=4)
				else:
					ax.fill_between(exon_x_region, exon_y_region_top, exon_y_region_bottom, facecolor=selected_color, alpha=1, edgecolor="none", zorder=4)

		displayed_trans_ID = re.sub('#Target_RBP_SDA|#Target_IDT|#Target_Ctrl|#Direct_RNA|#Target_SDA|_Sendai_03451_hICset4_D7','',trans_ID)
		if re.findall('ENST', displayed_trans_ID):
			displayed_trans_ID = displayed_trans_ID +' (%s)'%tag_dict[trans_ID]
		############# draw 5' and 3' end ###########
		[left_end, right_end] = decide_strand(strand_dict[trans_ID])
		if strand_dict[trans_ID] == '-':
			plt.text(-1*corrected_max_pos, line_space*(list_index+0.9)+exon_height, displayed_trans_ID, fontsize=font_size+1)
			plt.text(-1.03*corrected_max_pos, line_space*(list_index+1)-2, '5\'', fontsize=font_size+2)
		else:
			plt.text(0, line_space*(list_index+0.9)+exon_height, displayed_trans_ID, fontsize=font_size+1)
			plt.text(-0.03*corrected_max_pos, line_space*(list_index+1)-2, '5\'', fontsize=font_size+2)
		#plt.text(1.1*corrected_max_pos, line_space*(list_index+1)-3, right_end, fontsize=11)

#print ('corrected_max_pos', corrected_max_pos)

if strand_dict[trans_ID] == '-':
	ax.set_xlim(-1.1*corrected_max_pos, 0.1*corrected_max_pos)
	plt.text(-0.55*corrected_max_pos, line_space*(num_trans+1), gene_name, fontsize=font_size+2)
else:
	plt.text(0.45*corrected_max_pos, line_space*(num_trans+1), gene_name, fontsize=font_size+2)
	ax.set_xlim(-0.1*corrected_max_pos, 1.1*corrected_max_pos)
ax.set_ylim(0, line_space*(num_trans+1.5))
#plt.axis('equal')
plt.axis('off')
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
plt.margins(0,0)
plt.show()

file_dir = os.path.dirname(bed_list_name)
outf_name_pdf = f"{file_dir}/{gene_name}_{trans_ID_list[-1]}_isoform_structure_withORF_intron_shrinkage_{intron_shrinkage_fold}.pdf"
plt.savefig(outf_name_pdf,dpi=300, pad_inches=0)

outf_name_png = f"{file_dir}/{gene_name}_{trans_ID_list[-1]}_isoform_structure_withORF_intron_shrinkage_{intron_shrinkage_fold}.png"
plt.savefig(outf_name_png, dpi=300, pad_inches=0)
