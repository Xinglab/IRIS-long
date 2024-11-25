import os,re,sys
import time
from collections import defaultdict

netMHCpan_dir = sys.argv[1].rstrip('/')
protein_inf_name = sys.argv[2]
HLA_str = sys.argv[3]
window_size = sys.argv[4]
outf_dir = sys.argv[5].rstrip('/')
de_inf = sys.argv[6]
spe_inf = sys.argv[7]

#### generated differential isoform derived protein fasta #####
outf_pro = open("%s/6_1_Tumor_vs_normal_DE_prevalence_protein.fasta" % outf_dir,"w")
iso_dict = defaultdict()
with open(de_inf, "r") as DE_inf:
	for index,line in enumerate(DE_inf):
		if index == 0: continue
		arr = line.strip().split("\t")
		trans_ID = arr[0]
		iso_dict[trans_ID] = "Tumor-enriched"

with open(spe_inf, "r") as pre_inf:
	for index,line in enumerate(pre_inf):
		if index == 0: continue
		arr = line.strip().split("\t")
		trans_ID = arr[0]
		iso_dict[trans_ID] = "Tumor-specific"

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
HLA_list = HLA_str.split(",")
outf_sh = open("%s/6_1_perform_netMHCpan.sh" % outf_dir, "w")
outf_sh.write(f"#!/bin/bash\n#SBATCH -J netMHCpan\n#SBATCH -t 120:00:00\n#SBATCH -n 4\n#SBATCH --mem=120G\n#SBATCH --export=ALL\n#SBATCH --array=1-{len(HLA_list)}\n")
outf_sh.write("s=$(sed -n \"${SLURM_ARRAY_TASK_ID}p\" %s/6_1_perform_netMHCpan_cmd)\neval \"$s\"\n" % outf_dir)
#outf_sh.write("export s=`sed -n ${SLURM_ARRAY_TASK_ID}p %s/6_1_perform_netMHCpan_cmd`\necho $s\n$s\n" % outf_dir)
outf_sh.close()

outf_sh_array = open("%s/6_1_perform_netMHCpan_cmd" % outf_dir, "w")
outf_name_list = []
for each_HLA_str in HLA_list:
	outf_sh_array.write("%s/netMHCpan -inptype 0 -f %s/6_1_Tumor_vs_normal_DE_prevalence_protein.fasta -BA -t 2 -s -l %s -a %s > %s/6_1_TCR_temp_out_%s.txt\n" % (netMHCpan_dir, outf_dir, window_size, each_HLA_str, outf_dir, each_HLA_str))
	outf_name_list.append("%s/6_1_TCR_temp_out_%s.txt" % (outf_dir, each_HLA_str))
outf_sh_array.close()	

command_1 = f"sbatch {outf_dir}/6_1_perform_netMHCpan.sh"
print(command_1)
os.system(command_1)


##### check if files are completed or not #####
# Wait for job completion based on file size stability
job_completed = False
job_completed_dict = defaultdict(lambda: 0)
max_wait_time = 12 * 3600  # 12 hours
check_interval = 1800  # 30 minutes
waited_time = 0
previous_size_dict = defaultdict(lambda: -1)
current_size_dict = defaultdict()

while not job_completed and waited_time < max_wait_time:
	for each_file in outf_name_list:
		key_name = each_file.split("/")[-1].split(".")[0]
		if os.path.exists(each_file):
			current_size_dict[key_name] = os.path.getsize(each_file)
			if current_size_dict[key_name] == previous_size_dict[key_name]:  # File size hasn't changed
				job_completed_dict[key_name] = 1
			else:
				previous_size_dict[key_name] = current_size_dict[key_name]  # Update previous size
		else:
			previous_size_dict[key_name] = -1  # Reset if file doesn't exist

	# Check if all jobs are completed
	print(job_completed_dict)
	if sum(job_completed_dict.values()) == len(outf_name_list):
		job_completed = True
	else:
		time.sleep(check_interval)
		waited_time += check_interval

if not job_completed:
	print("Job did not complete within the expected time or file size stopped updating.")
	sys.exit(1)


###### Once all files are generated #####
os.system("cat %s/6_1_TCR_temp_out_*.txt > %s/6_1_TCR_temp_out.txt" % (outf_dir, outf_dir))

############## process netMHCpan result #############
peptide_info_dict = defaultdict(lambda: [])
info_dict = defaultdict(lambda: '')
with open('%s/6_1_TCR_temp_out.txt' % outf_dir, 'r') as inf_1:
	for index, line in enumerate(inf_1):
		arr = line.strip().split()
		if len(arr) == 0: continue
		if arr[-1] in ["SB","WB"]:
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

