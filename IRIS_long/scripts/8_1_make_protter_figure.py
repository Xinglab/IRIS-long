import os,re,sys
from collections import defaultdict
import requests

in_trans_ID = sys.argv[1]
in_gene_name = sys.argv[2]
protein_inf = sys.argv[3]
outf_dir = sys.argv[4]
high_specificity_cutoff = float(sys.argv[5])
low_specificity_cutoff = int(1)

server = "http://wlab.ethz.ch/protter/create?"

trans2pro_seq_dict = defaultdict()
trans_ID = ''
with open(protein_inf, 'r') as protein_inf:
	for line in protein_inf:
		if line.startswith('>'):
			trans_ID = line.strip().split('|')[2].split('_')[0]
		else:
			trans2pro_seq_dict[trans_ID] = line.strip()

final_url_list = [server]
pos2score_dict = defaultdict(lambda: 0)

def convert_to_str(pos_list):
	if len(pos_list) == 0: return ''
	region_list = []
	prev_region = pos_list[0]+'-'+pos_list[0]
	for i in range(1,len(pos_list)):
		if int(pos_list[i]) == int(prev_region.split('-')[1])+1:
			prev_region = prev_region.split('-')[0]+'-'+pos_list[i]
		else:
			region_list.append(prev_region)
			prev_region = pos_list[i]+'-'+pos_list[i]
	region_list.append(prev_region)
	return ','.join(region_list)


with open("%s/5_5_Summarized_CAR_T_prioritized_targets_final.txt" % outf_dir,'r') as inf:
	for index,line in enumerate(inf):
		arr = line.strip().split('\t')
		if index == 0: continue
		if arr[0] != in_trans_ID: continue
		### sequence ###
		if arr[16] != '-':
			final_url_list.append("up="+arr[16])
		else:
			final_url_list.append("seq="+trans2pro_seq_dict[arr[0]])
		### N terminal ###
		if arr[22] in ['Annotated_cell_surface:TM', 'Identified by both: consistent']:
			if re.findall("TMD",arr[20]):
				if arr[19].split(':')[1].startswith('1-'):
					final_url_list.append("nterm=extra")
				else:
					final_url_list.append("nterm=intra")
			elif re.findall('0-1', arr[20]):
				final_url_list.append("nterm=extra")
			else:
				final_url_list.append("nterm=intra")
			### TM pos ###
			final_url_list.append("tm="+arr[17].split(':')[1])
		else:
			final_url_list.append("nterm=extra")
		### get tumor specificity score each pos ###
		for i in range(len(arr[26].split(";"))):
			score = float(arr[27].split(";")[i])
			region = arr[26].split(";")[i]
			for i_pos in range(int(region.split('-')[0]), int(region.split('-')[1])+1):
				if score > pos2score_dict[str(i_pos)]:
					pos2score_dict[str(i_pos)] = score
		### target regions ###
		low_score_pos_list = []
		high_score_pos_list = []
		for each_target_region in arr[3].split(";"):
			for each_pos in range(int(each_target_region.split("-")[0]), int(each_target_region.split("-")[1])+1):
				if (pos2score_dict[str(each_pos)] > low_specificity_cutoff) and (pos2score_dict[str(each_pos)] <= high_specificity_cutoff):
					low_score_pos_list.append(str(each_pos))
				elif pos2score_dict[str(each_pos)] > high_specificity_cutoff:
					high_score_pos_list.append(str(each_pos))
		low_score_pos_str = convert_to_str(low_score_pos_list)
		high_score_pos_str = convert_to_str(high_score_pos_list)
		print(low_score_pos_str, high_score_pos_str)
		final_url_list.append(f"n:Target regions with tumor specificity score from {low_specificity_cutoff}-{high_specificity_cutoff},s:circ,fc:black,bc:orange=" + low_score_pos_str)
		final_url_list.append(f"n:Target regions with tumor specificity score greater than {high_specificity_cutoff},s:circ,fc:black,bc:red=" + high_score_pos_str)
		### parameters ###
		final_url_list.append("mc=lightsalmon&lc=blue&tml=numcount&numbers&legend")
		final_url = '&'.join(final_url_list)
		break
print (final_url)

outf_name = "%s/8_1_Protter_%s_%s" % (outf_dir, in_gene_name,in_trans_ID)
#with open(outf_name+'.png', 'wb') as f_png:
#	response = requests.get(final_url+'&format=png')
#	f_png.write(response.content)
with open(outf_name+'.pdf', 'wb') as f_pdf:
	response = requests.get(final_url+'&format=pdf')
	f_pdf.write(response.content)
with open(outf_name+'.svg', 'wb') as f_svg:
	response = requests.get(final_url+'&format=svg')
	f_svg.write(response.content)

