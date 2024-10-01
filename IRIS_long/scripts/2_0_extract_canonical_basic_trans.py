import os,re,sys
from collections import defaultdict

gtf_inf_path = sys.argv[1]

file_dir = os.path.dirname(os.path.realpath(__file__))
version = gtf_inf_path.split("/")[-1].split(".")[1]
basic_outf_path = os.path.dirname(gtf_inf_path)+f"/gencode.{version}.basic_transcript.txt"
canonical_outf_path = os.path.dirname(gtf_inf_path)+f"/gencode.{version}.canonical_transcript.txt"

def get_right_ID(raw_ID):
	if re.findall("_PAR_Y", raw_ID):
		return raw_ID.split(".")[0]+"-PAR-Y"
	else:
		return raw_ID.split(".")[0]

GeneID2name_dict = defaultdict()
GeneID2cano_dict = defaultdict()
GeneID2basic_dict = defaultdict(lambda: [])
with open(gtf_inf_path,"r") as gtf_inf:
	for line in gtf_inf:
		if line.startswith('#'): continue
		arr = line.strip().split("\t")
		if arr[2] != "transcript": continue
		trans_ID = get_right_ID(re.findall(r'transcript_id "(.+?)";',arr[8])[0])
		gene_ID = get_right_ID(re.findall(r'gene_id "(.+?)";',arr[8])[0])
		gene_name = re.findall(r'gene_name "(.+?)";',arr[8])[0]
		tag_list = re.findall(r'tag "(.+?)";',arr[8])
		GeneID2name_dict[gene_ID] = gene_name
		if "Ensembl_canonical" in tag_list:
			GeneID2cano_dict[gene_ID] = trans_ID
		if "basic" in tag_list: 
			GeneID2basic_dict[gene_ID].append(trans_ID)

canonical_outf = open(canonical_outf_path, "w")
canonical_outf.write("Gene_ID\tGeneSymbol\tTranscript_ID\n")
for each_gene in GeneID2cano_dict:
	canonical_outf.write(f"{each_gene}\t{GeneID2name_dict[each_gene]}\t{GeneID2cano_dict[each_gene]}\n")
canonical_outf.close()

basic_outf = open(basic_outf_path, "w")
basic_outf.write("Gene_ID\tGeneSymbol\tTranscript_ID\n")
for each_gene in GeneID2basic_dict:
	basic_str = ";".join(GeneID2basic_dict[each_gene])
	basic_outf.write(f"{each_gene}\t{GeneID2name_dict[each_gene]}\t{basic_str}\n")
basic_outf.close()

