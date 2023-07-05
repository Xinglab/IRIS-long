import os,re,sys
from collections import defaultdict

netMHCpan_dir = sys.argv[1].rstrip('/')
protein_inf_name = sys.argv[2]
HLA_str = sys.argv[3]
window_size = sys.argv[4]
outf_dir = sys.argv[5].rstrip('/')

#### generated differential isoform derived protein fasta #####
outf_pro = open("%s/6_1_Tumor_vs_normal_DE_prevalence_protein.fasta" % outf_dir,"w")
iso_dict = defaultdict()
with open("%s/../3_1_Tumor_vs_normal_DE_test.txt" % outf_dir,"r") as DE_inf:
	for index,line in enumerate(DE_inf):
		if index == 0: continue
		arr = line.strip().split("\t")
		trans_ID = arr[0]
		iso_dict[trans_ID] = 1

with open("%s/../3_2_Tumor_vs_normal_prevalence.txt" % outf_dir,"r") as pre_inf:
	for index,line in enumerate(pre_inf):
		if index == 0: continue
		arr = line.strip().split("\t")
		trans_ID = arr[0]
		iso_dict[trans_ID] = 1

flag = 1
with open(protein_inf_name, "r") as pro_inf:
	for line in pro_inf:
		if line.startswith(">"):
			flag = 0
			trans_ID = line.strip().split("|")[2].split("_")[0]
			if trans_ID in iso_dict:
				outf_pro.write(line)
				flag = 1
		else:
			if flag == 1:
				outf_pro.write(line)
outf_pro.close()

############## generate netMHCpan bash file ##########
outf_sh = open("%s/6_1_perform_netMHCpan.sh" % outf_dir, "w")
outf_sh.write("#!/bin/bash\n#SBATCH -J netMHCpan\n#SBATCH -o netMHCpan.out\n#SBATCH -e netMHCpan.err\n#SBATCH -t 120:00:00\n#SBATCH -n 8\n#SBATCH --mem=40G\n#SBATCH --export=ALL\n")
outf_sh.write("%s/netMHCpan -inptype 0 -f %s/6_1_Tumor_vs_normal_DE_prevalence_protein.fasta -BA -t 2 -s -l %s -a %s > %s/6_1_TCR_temp_out.txt\n" % (netMHCpan_dir, outf_dir, window_size, HLA_str, outf_dir))
outf_sh.close()	
command_1 = f"bash {outf_dir}/6_1_perform_netMHCpan.sh"
print(command_1)
os.system(command_1)


############## process netMHCpan result #############
peptide_info_dict = defaultdict(lambda: [])
info_dict = defaultdict(lambda: '')
with open('%s/6_1_TCR_temp_out.txt' % outf_dir, 'r') as inf_1:
	for index, line in enumerate(inf_1):
		arr = line.strip().split()
		if len(arr) == 0: continue
		if arr[-1] in ["SB"]:
			[pos, HLA_type, peptide, protein, rank_EL, rank_BA, Aff, level] = [arr[0], arr[1], arr[2], arr[10].split('_')[1], arr[12], arr[14], arr[15], arr[-1]]
			key = peptide+'_'+HLA_type
			value = "\t".join([peptide, HLA_type, rank_EL, rank_BA, Aff, level])
			info_dict[key] = value
			value2 = ",".join([HLA_type, rank_EL, rank_BA, Aff, level])
			if value2 not in peptide_info_dict[peptide]:
				peptide_info_dict[peptide].append(value2)

outf = open("%s/6_1_Summarized_TCR_netMHCpan_out.txt" % outf_dir, "w")
outf.write("Peptide\tHLA_type\t%Rank_EL\t%Rank_BA\tAff_nM\tLevel\n")
for each_key in info_dict:
	outf.write(info_dict[each_key]+'\n')
outf.close()

outf2 = open("%s/6_1_Summarized_TCR_netMHCpan_reshaped_out.txt" % outf_dir, "w")
outf2.write("Peptide\tHLA_type\t%Rank_EL\t%Rank_BA\tAff_nM\tLevel\n")
for each_key in peptide_info_dict:
	a_list = [each_value.split(',') for each_value in peptide_info_dict[each_key]]
	outf2.write(each_key)
	for each_tuple in zip(*a_list):
		outf2.write('\t'+';'.join(list(each_tuple)))
	outf2.write('\n')
outf2.close()

