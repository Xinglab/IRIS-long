import os,re,sys
from collections import defaultdict

dir_path = os.path.dirname(os.path.realpath(__file__))
inf_name = sys.argv[1] #'CD44_sample_list_N2_R0_updated.gtf'
reference_gtf = sys.argv[3]

if re.findall("encode", inf_name.split('/')[-1]):
	key_word = 'gencode'
else:
	key_word = inf_name.split('/')[-1].split('_')[0]

if len(sys.argv) > 4:
	outf_path = sys.argv[4]
else:
	outf_path = '/'.join(os.path.abspath(inf_name).split('/')[0:-1])

if len(sys.argv) > 2:
	inf_abundance_name = sys.argv[2] # 'abundance.esp'
else:
	inf_abundance_name = 'NA'

outf_name = outf_path+'/samples_BedGraph.bed'
outf = open(outf_name, 'w')
outf.write('track name="%s" visibility=dense\n' % key_word)

def get_right_ID(raw_ID):
	if re.findall('_PAR_Y', raw_ID):
		return_ID = raw_ID.split('.')[0]+'-PAR-Y'
	else:
		return_ID = raw_ID.split('.')[0]
	return return_ID

def convert_list_to_str(trans_range,exon_list,trans_anno):
	chrom = trans_range.split('$')[0]
	strand = trans_range.split('$')[1]
	chromStart = int(trans_range.split('$')[2]) - 1 # 0-based
	chromEnd = int(trans_range.split('$')[3])
	thickStart = chromStart
	thickEnd = chromStart
	itemRgb = '0'
	name = exon_list[0].split('$')[0] #ENST_ID
	score = 0
	blockCount = len(exon_list)
	blockSizes = []
	blockStarts = []
	if len(exon_list) > 1:
		first_exon_start = int(exon_list[0].split('$')[3]) - 1
		second_exon_start = int(exon_list[1].split('$')[3]) - 1
		if first_exon_start > second_exon_start:
			exon_list = exon_list[::-1]
	for each_exon in exon_list:
		exon_start = int(each_exon.split('$')[3]) - 1
		exon_end = int(each_exon.split('$')[4])
		blockSizes.append(str(exon_end - exon_start))
		blockStarts.append(str(exon_start - chromStart))
	res = chrom+'\t'+str(chromStart)+'\t'+str(chromEnd)+'\t'+name+'\t'+str(score)+'\t'+strand+'\t'+str(thickStart)+'\t'+str(thickEnd)+'\t'+itemRgb+'\t'+str(blockCount)+'\t'+','.join(blockSizes)+'\t'+','.join(blockStarts)
	return res

#### load some annotation ####
trans2gene_dict = defaultdict()
with open(reference_gtf, "r") as gtf_inf:
	for line in gtf_inf:
		if line.startswith('#'): continue
		arr = line.strip().split('\t')
		if arr[2] != 'transcript': continue
		ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\"',arr[8])[0])
		ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\"',arr[8])[0])
		if ENST_ID not in trans2gene_dict:
			trans2gene_dict[ENST_ID] = ENSG_ID

if inf_abundance_name != 'NA':
	with open(inf_abundance_name, 'r') as abun_inf:
		for index, line in enumerate(abun_inf):
			if index == 0: continue
			arr = line.strip().split('\t')
			ENST_ID = get_right_ID(arr[0])
			ENSG_ID = get_right_ID(arr[2])
			if ENST_ID not in trans2gene_dict:
				trans2gene_dict[ENST_ID] = ENSG_ID

######################
key_list_dict = defaultdict(lambda: [])
transcript_annotate_dict = defaultdict()
transcript_range_dict = defaultdict()
trans_list = []

#file_len = int(os.popen('wc -l %s' % inf_name).readlines()[0].split(' ')[0])
inf = open(inf_name,'r')
for i,line in enumerate(inf):
	if line.startswith('#'): continue
	arr = line.strip().split('\t')
	if arr[2] == 'transcript':
		ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\"',arr[8])[0])
		trans_list.append(ENST_ID)
		transcript_annotate_dict[ENST_ID] = arr[1]
		transcript_range_dict[ENST_ID] = arr[0]+'$'+arr[6]+'$'+str(int(arr[3]))+'$'+str(arr[4])
	elif arr[2] == 'exon':
		ENST_ID = get_right_ID(re.findall('transcript_id \"(.+?)\"',arr[8])[0])
		if re.findall("gene_id", arr[8]):
			ENSG_ID = get_right_ID(re.findall('gene_id \"(.+?)\"',arr[8])[0])
		elif ENST_ID in trans2gene_dict:
			ENSG_ID = trans2gene_dict[ENST_ID]
		else:
			continue
		key = ENST_ID+'|'+ENSG_ID+'$'+arr[0]+'$'+arr[6]+'$'+str(int(arr[3]))+'$'+str(arr[4])
		key_list_dict[ENST_ID].append(key)
inf.close()

for each_trans_ID in trans_list:
	out_str = convert_list_to_str(transcript_range_dict[each_trans_ID],key_list_dict[each_trans_ID],transcript_annotate_dict[each_trans_ID])
	outf.write(out_str+'\n')
outf.close()	
