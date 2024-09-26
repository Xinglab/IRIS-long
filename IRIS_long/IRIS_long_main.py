import os,re,sys
from collections import defaultdict
import configargparse
from statsmodels.stats.multitest import multipletests

__author__ = 'Yang Xu'
#__version__ = config.CURRENT_VERSION
__email__ = 'yangax@pennmedicine.upenn.edu'

def parse_args():
	parser = configargparse.ArgParser(description='IRIS-long workflow')
	subparsers = parser.add_subparsers(title='Sub-commands', description='Combine\nPreprocess\nFigure\nDiffTest\nTranslation\nCAR_T\nTCR\nSpecificity\n', dest='subcommand')
	
	# create the parser for the 'Combine' command
	parser_combine = subparsers.add_parser('Combine', help='combine ESPRESSO outputs from different batch')
	parser_combine.add_argument('-l', '--gtf_list', dest='gtf_list', type=str, help='list of GTF files that needs to be combined', required=True)
	parser_combine.add_argument('-ad', '--allowed_dist', dest='allowed_dist', type=int, help='allowed distance (bp) for each end to determine two novel transcripts is the same one', default = 50)
	#parser_combine.add_argument('-oe', '--output_abundance', dest='output_abundance', type=str, help='filename of the output combined abundance matrix', required=True)
	#parser_combine.add_argument('-og', '--output_gtf', dest='output_gtf', type=str, help='filename of the output combined gtf file', required=True)
	parser_combine.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'Preprocess' command
	parser_preprocess = subparsers.add_parser('Preprocess', help='generate CPM file as well as isoform proportion matrix')    
	parser_preprocess.add_argument('-ig', '--input_gtf', dest='input_gtf', type=str, help='Input gtf file', required=True)
	parser_preprocess.add_argument('-ia', '--input_abundance', dest='input_abundance', type=str, help='Input abundance file', required=True)
	parser_preprocess.add_argument('-rg', '--ref_gtf', dest='ref_gtf', type=str, help='Reference Gencode gtf', default = './')
	parser_preprocess.add_argument('-nm', '--normalized_mode', dest='normalized_mode', type=str, choices=['SAM','SELF'], help="Choose normalization mode from ['SAM','SELF']", required=True)
	parser_preprocess.add_argument('-fs', '--folder_sam', dest='folder_sam', type=str, help='directory of corresponding sam files, sample names need to match abundance matrix', default = './')
	parser_preprocess.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'Figure' command
	parser_figure = subparsers.add_parser('Figure', help='generate barplot and isoform structure plot')
	parser_figure.add_argument('-ip', '--isoform_proportion_inf', dest='isoform_proportion_inf', type=str, help='Isoform proportion infile, e.g. samples_abundance_combined_CPM_proportion.txt', required=True)
	parser_figure.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM file, e.g. samples_abundance_combined_CPM.txt', required=True)
	parser_figure.add_argument('-gi', '--group_info_inf', dest='group_info_inf', type=str, help='Group information file', required=True)
	parser_figure.add_argument('-rt', '--required_trans_inf', dest='required_trans_inf', type=str, help='Required transcripts information file', required=True)
	parser_figure.add_argument('-be', '--bedgraph', dest='bedgraph', type=str, help='Generated bedgraph file for sample, e.g. samples_BedGraph.bed', required=True)
	parser_figure.add_argument('-gv', '--genome_version', dest='genome_version', type=str, choices=['GRCh38','GRCh37','hg38','hg19'], help="choose from ['GRCh38','GRCh37','hg38','hg19']", default='hg38')
	parser_figure.add_argument('-fi', '--figures', dest='figures', nargs='+', type=str, choices=['Isoform','Single_isoform','Structure'], default=['Isoform','Single_isoform','Structure'],  help='Choose from [Isoform, Single_isoform, Structure]', required=True)
	parser_figure.add_argument('-ci', '--CDS_inf', dest='CDS_inf', type=str, help="generated CDS file, e.g. 4_4_*_detailed_match_ID.txt", required=True)
	parser_figure.add_argument('-is', '--intron_shrinkage', dest='intron_shrinkage', type=int, help="Intron shrinkage fold in isoform structure figure", default=10)
	parser_figure.add_argument('-of', '--order', dest='order', type=str, help='order samples by isoform proportion', default='no')
	parser_figure.add_argument('-rg', '--ref_gtf', dest='ref_gtf', type=str, help='Reference Gencode gtf', default = './')
	parser_figure.add_argument('-ig', '--input_gtf', dest='input_gtf', type=str, help='Input gtf file', required=True)
	parser_figure.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'DiffTest' command
	parser_difftest = subparsers.add_parser('DiffTest', help='Differential tests between tumor and normal tissues')
	parser_difftest.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM file, e.g. samples_abundance_combined_CPM.txt', required=True)
	parser_difftest.add_argument('-tn', '--tumor_num', dest='tumor_num', type=int, help='Number of tumor samples', required=True)
	parser_difftest.add_argument('-ep', '--enriched_test_p', dest='enriched_test_p', type=float, default=0.05, help='Cutoff of p-value in DE test (default = 0.05)')
	parser_difftest.add_argument('-etc', '--enriched_test_tumor_cpm', dest='enriched_test_tumor_cpm', type=float, help='Cutoff of median CPM value in tumor samples, used to filter out lowly-expressed isoform (default = 5)', default = 5)
	parser_difftest.add_argument('-ef', '--enriched_test_fc', dest='enriched_test_fc', type=float, help='Cutoff of fold change between tumor and normal, used to decide DE isoform (default = 2)', default = 2)
	#parser_difftest.add_argument('-sp', '--specificity_test_p', dest='specificity_test_p', type=float, default=1e-6, help='Cutoff of p-value in prevalence test (default = 1e-6)')
	parser_difftest.add_argument('-stc', '--specificity_test_tumor_cpm', dest='specificity_test_tumor_cpm', type=float, help='Cutoff of CPM value in tumor samples, used to decide whether an isoform is considered as expressed (default = 5)', default = 5)
	parser_difftest.add_argument('-snc', '--specificity_test_tissue_cpm', dest='specificity_test_tissue_cpm', type=float, help='Cutoff of CPM value in tissue samples, used to decide whether an isoform is considered as expressed (default = 1)', default = 1)
	parser_difftest.add_argument('-stp', '--specificity_test_tumor_percentage', dest='specificity_test_tumor_percentage', type=float, help='Minimum precentage of tumor samples that express given transcript (default = 50%, which is 0.5)', default = 0.5)
	parser_difftest.add_argument('-snp', '--specificity_test_tissue_percentage', dest='specificity_test_tissue_percentage', type=float, help='Maximum precentage of tissue samples that express given transcript (default = 10%, which is 0.1)', default = 0.1)
	parser_difftest.add_argument('-sen', '--specificity_test_exclude_tissue', dest='specificity_test_exclude_tissue', type=str, help='Excluded tissues in tumor-specificity transcript identification, separated by comma', default = 'Testis')
	parser_difftest.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'Translation' command
	parser_translate = subparsers.add_parser('Translation', help='translate transcripts into protein sequence')
	parser_translate.add_argument('-mo', '--mode', dest='mode', type=str, default = 'long-read', help='Long-read RNA-seq data mode')
	parser_translate.add_argument('-tg', '--trans_gtf', dest='trans_gtf', type=str, help='generated gtf file, e.g. samples_updated_combined.gtf', required=True)
	parser_translate.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM file, e.g. samples_abundance_combined_CPM.txt', required=True)
	parser_translate.add_argument('-gv', '--genome_version', dest='genome_version', type=str, choices=['GRCh38','GRCh37','hg38','hg19'], help="choose from ['GRCh38','GRCh37','hg38','hg19']", required=True)
	parser_translate.add_argument('-rg', '--ref_gtf', dest='ref_gtf', type=str, help='reference gencode annotation, e.g. gencode.v39.annotation.gtf', default = './')
	parser_translate.add_argument('-of', '--out_file', dest='out_file', type=str, help='prefix of the name of output file', required=True)
	parser_translate.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'CAR-T prediction' command
	parser_car_t = subparsers.add_parser('CAR_T', help='CAR_T target prediction')
	parser_car_t.add_argument('-td', '--tmhmm_dir', dest='tmhmm_dir', type=str, help='file path of TMHMM tool (directory is needed)', required=True)
	parser_car_t.add_argument('-tn', '--tumor_num', dest='tumor_num', type=int, help='Number of tumor samples', required=True)
	parser_car_t.add_argument('-pi', '--protein_inf', dest='protein_inf', type=str, help='Generated protein fasta file', required=True)
	parser_car_t.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM file, e.g. samples_abundance_combined_CPM.txt', required=True)
	parser_car_t.add_argument('-ip', '--isoform_proportion_inf', dest='isoform_proportion_inf', type=str, help='Isoform proportion file, e.g. samples_abundance_combined_CPM_proportion.txt', required=True)
	parser_car_t.add_argument('-gv', '--genome_version', dest='genome_version', type=str, choices=['GRCh38','GRCh37','hg38','hg19'], help="choose from ['GRCh38','GRCh37','hg38','hg19']", required=True)
	parser_car_t.add_argument('-ss', '--specificity_score', dest='specificity_score', type=float, help='cutoff of specificity_score (default = 1)', default = 1)
	parser_car_t.add_argument('-tc', '--tissue_cpm', dest='tissue_cpm', type=float, help='cutoff of (maximum tolerable) CPM of transcripts encode given peptide in tissue samples (default = 10)', default = 10)
	parser_car_t.add_argument('-tcp', '--tissue_percentage', dest='tissue_percentage', type=float, help='maximum tolerable percentage of tissues that are allowed to have transcript higher than the given CPM expression threshold (default = 20%, which is 0.2)', default = 0.2)
	parser_car_t.add_argument('-aic', '--annotated_isoform_contri_inf', dest='annotated_isoform_contri_inf', type=str, help='file generated before, which ends with \"_annotated_isoform_contribution.txt\"', required=True)
	parser_car_t.add_argument('-tci', '--trans_CDS_inf', dest='trans_CDS_inf', type=str, help='file generated before, format of which is like \"4_4_*_detailed_match_ID.txt\"', required=True)
	#parser_car_t.add_argument('-oic', '--other_isoform_contri_inf', dest='other_isoform_contri_inf', type=str, help='file generated before, which ends with \"_proportion_only_focus_others.txt\"', required=True)
	parser_car_t.add_argument('-rg', '--ref_gtf', dest='ref_gtf', type=str, help='reference gencode annotation, e.g. gencode.v39.annotation.gtf', default = './')
	parser_car_t.add_argument('-gf', '--gencode_fasta', dest='gencode_fasta', type=str, help='reference gencode translation.fasta, e.g. gencode.v39.pc_translations.fa', default = './')
	parser_car_t.add_argument('-dt', '--de_trans_inf', dest='de_trans_inf', type=str, help='Tumor-enriched transcripts, e.g. 3_1_Tumor_vs_normal_DE_test.txt', required=True)
	parser_car_t.add_argument('-st', '--spe_trans_inf', dest='spe_trans_inf', type=str, help='Tumor-specific transcripts, e.g. 3_2_Tumor_vs_normal_prevalence.txt', required=True)
	parser_car_t.add_argument('-ws', '--window_size', dest='window_size', type=int, help='Window size (default = 9 AAs)', default = 9)
	parser_car_t.add_argument('-of', '--out_file', dest='out_file', type=str, help='prefix of the name of output file', required=True)
	parser_car_t.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'TCR prediction' command
	parser_tcr = subparsers.add_parser('TCR', help='TCR target prediction')
	parser_tcr.add_argument('-nd', '--netMHCpan_dir', dest='netMHCpan_dir', type=str, help='file path of netMHCpan_dir tool (directory is needed)', required=True)
	parser_tcr.add_argument('-hs', '--HLA_str_inf', dest='HLA_str_inf', type=str, help='File containing HLA alleles information, first column is sample, and second column is interesed HLA allele that separated by comma, e.g. HLA-A02:01,HLA-A01:01', required=True)
	parser_tcr.add_argument('-tn', '--tumor_num', dest='tumor_num', type=int, help='Number of tumor samples', required=True)
	parser_tcr.add_argument('-pi', '--protein_inf', dest='protein_inf', type=str, help='Generated protein fasta file, such as 4_4_XXX_PC.fasta', required=True)
	parser_tcr.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM file, e.g. samples_abundance_combined_CPM.txt', required=True)
	parser_tcr.add_argument('-ip', '--isoform_proportion_inf', dest='isoform_proportion_inf', type=str, help='Isoform proportion file, e.g. samples_abundance_combined_CPM_proportion.txt', required=True)
	parser_tcr.add_argument('-gv', '--genome_version', dest='genome_version', type=str, choices=['GRCh38','GRCh37','hg38','hg19'], help="choose from ['GRCh38','GRCh37','hg38','hg19']", required=True)
	parser_tcr.add_argument('-ss', '--specificity_score', dest='specificity_score', type=float, help='cutoff of specificity_score (default = 3)', default = 3)
	parser_tcr.add_argument('-tc', '--tissue_cpm', dest='tissue_cpm', type=float, help='cutoff of (maximum tolerable) CPM of transcripts encode given peptide in tissue samples (default = 10)', default = 10)
	parser_tcr.add_argument('-tcp', '--tissue_percentage', dest='tissue_percentage', type=float, help='maximum tolerable percentage of tissues that are allowed to have peptide higher than the given CPM expression threshold (default = 20%, which is 0.2)', default = 0.2)
	parser_tcr.add_argument('-ba', '--binding_affi', dest='binding_affi', type=float, help='cutoff of binding affinity between HLA complex and peptide (default = 500)', default = 500)
	parser_tcr.add_argument('-aic', '--annotated_isoform_contri_inf', dest='annotated_isoform_contri_inf', type=str, help='file generated before, which is like \"_annotated_isoform_contribution.txt\"', required=True)
	parser_tcr.add_argument('-tci', '--trans_CDS_inf', dest='trans_CDS_inf', type=str, help='file generated before, which is like  \"4_4_*_detailed_match_ID.txt\"', required=True)
	parser_tcr.add_argument('-rg', '--ref_gtf', dest='ref_gtf', type=str, help='reference gencode annotation, e.g. gencode.v39.annotation.gtf', default = './')
	parser_tcr.add_argument('-dt', '--de_trans_inf', dest='de_trans_inf', type=str, help='Tumor-enriched transcripts, e.g. 3_1_Tumor_vs_normal_DE_test.txt', required=True)
	parser_tcr.add_argument('-st', '--spe_trans_inf', dest='spe_trans_inf', type=str, help='Tumor-specific transcripts, e.g. 3_2_Tumor_vs_normal_prevalence.txt', required=True)
	parser_tcr.add_argument('-ws', '--window_size', dest='window_size', type=int, help='Window size (default = 9 AAs)', default = 9)
	parser_tcr.add_argument('-of', '--out_file', dest='out_file', type=str, help='prefix of the name of output file', required=True)
	parser_tcr.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'Specificity' command
	parser_specificity = subparsers.add_parser('Specificity', help='Sliding window to find peptide with the highest tumor-specificity')
	parser_specificity.add_argument('-ti', '--transcript_ID', dest='transcript_ID', type=str, help='Interested transcript ID', required=True)
	parser_specificity.add_argument('-tn', '--tumor_num', dest='tumor_num', type=int, help='Number of tumor samples', required=True)
	parser_specificity.add_argument('-pi', '--protein_inf', dest='protein_inf', type=str, help='Generated protein fasta file, which is like \"4_4_*_PC.fasta\"', required=True)
	parser_specificity.add_argument('-ic', '--isoform_cpm_inf', dest='isoform_cpm_inf', type=str, help='Isoform CPM file, e.g. samples_abundance_combined_CPM.txt', required=True)
	parser_specificity.add_argument('-cs', '--cell_surface_inf', dest='cell_surface_inf', type=str, help='Generated cell surface proteins file, which is like \"5_3_*_high_confidence.txt\"', required=True)
	parser_specificity.add_argument('-ws', '--window_size', dest='window_size', type=int, help='Window size (default = 9 AAs)', default = 9)
	parser_specificity.add_argument('-ss', '--start_site', dest='start_site', type=int, help='Starting position of visualized protein region (default shows the 100 AAs region with the highest tumor-specificity)', default = 0)
	parser_specificity.add_argument('-es', '--end_site', dest='end_site', type=int, help='Ending position of visualized protein region (default shows the 100 AAs region with the highest tumor-specificity)', default = 999999)
	parser_specificity.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	# create the parser for the 'Topology' command
	parser_topology = subparsers.add_parser('Topology', help='Making protein topology figure using Protter tool')
	parser_topology.add_argument('-ti', '--transcript_ID', dest='transcript_ID', type=str, help='Interested transcript ID', required=True)
	parser_topology.add_argument('-gs', '--gene_symbol', dest='gene_symbol', type=str, help='Interested gene name', required=True)
	parser_topology.add_argument('-pi', '--protein_inf', dest='protein_inf', type=str, help='Generated protein fasta file, which is like \"4_4_*_PC.fasta\"', required=True)
	parser_topology.add_argument('-sc', '--score_cutoff', dest='score_cutoff', type=float, help='Specificity score cutoff (default = 3)', default = 3)
	parser_topology.add_argument('-od', '--outf_dir', dest='outf_dir', type=str, help='output directory', required=True)

	args = parser.parse_args()
	return parser, args


def main():
	parser, args = parse_args()
	#outf_dir = args.outf_dir
	dir_path = os.path.dirname(os.path.realpath(__file__))
	subcommand = args.subcommand
	if subcommand == "Combine":
		Sub_combine(dir_path, args)
	elif subcommand == "Preprocess":
		Sub_process(dir_path, args)
	elif subcommand == 'Figure':
		Sub_figure(dir_path, args)
	elif subcommand == 'DiffTest':
		Sub_difftest(dir_path, args)
	elif subcommand == 'Translation':
		Sub_translation(dir_path, args)
	elif subcommand == 'CAR_T':
		Sub_CAR_T(dir_path, args)
	elif subcommand == 'TCR':
		Sub_TCR(dir_path, args)
	elif subcommand == 'Specificity':
		Sub_specificity(dir_path, args)
	elif subcommand == 'Topology':
		Sub_topology(dir_path, args)


def Sub_combine(dir_path, args):
	gtf_list = args.gtf_list
	outf_dir = args.outf_dir.rstrip('/')
	allowed_dist = args.allowed_dist
	# 0.1 match novel ESPRESSO isoforms
	cmd_combine_1 = f"python {dir_path}/scripts/0_1_match_isoform_vs_novel_gtf.py {gtf_list} {allowed_dist} {outf_dir}/1_match_isoform.txt {outf_dir}/1_all_isoform.txt"
	print (cmd_combine_1)
	os.system(cmd_combine_1)
	# 0.2 merge abundance esp file
	cmd_combine_2 = f"python {dir_path}/scripts/0_2_merge_exp_matrix.py {outf_dir}/1_match_isoform.txt {outf_dir}/1_all_isoform.txt {gtf_list} {outf_dir}/samples_abundance_combined.txt"
	print (cmd_combine_2)
	os.system(cmd_combine_2)
	# 0.3 merge gtf file
	cmd_combine_3 = f"python {dir_path}/scripts/0_3_merge_gtf.py {outf_dir}/1_match_isoform.txt {outf_dir}/1_all_isoform.txt {gtf_list} {outf_dir}/samples_updated_combined.gtf"
	print (cmd_combine_3)
	os.system(cmd_combine_3)


def Sub_process(dir_path, args):
	input_gtf = args.input_gtf
	input_abundance = args.input_abundance
	ref_gtf = args.ref_gtf
	folder_sam = args.folder_sam
	normalized_mode = args.normalized_mode
	outf_dir = args.outf_dir.rstrip('/')
	# 1.1 convert gtf to bed file
	cmd_process_1 = f"python {dir_path}/scripts/1_1_gtf2bed_1_to_1_based.py {input_gtf} {input_abundance} {ref_gtf} {outf_dir}"
	print (cmd_process_1)
	os.system(cmd_process_1)
	# 1.2 calculate CPM based on sam/ban files
	cmd_process_2 = f"python {dir_path}/scripts/1_2_calculate_CPM_based_on_bam.py {input_abundance} {folder_sam} {outf_dir} {normalized_mode}"
	print (cmd_process_2)
	os.system(cmd_process_2)
	# 1.3 merge to gene level, and generate files to filter out bad-quality genes
	out_file_name = re.sub(".txt|.esp", f"_CPM.txt", input_abundance.split('/')[-1])
	cmd_process_3 = f"python {dir_path}/scripts/1_3_merge_to_gene_CPM.py {outf_dir}/{out_file_name}"
	print (cmd_process_3)
	os.system(cmd_process_3)
	# 1.4 calculate isoform proportion
	cmd_process_4 = f"python {dir_path}/scripts/1_4_calculate_isoform_proportion_vector.py {outf_dir}/{out_file_name}"
	print (cmd_process_4)
	os.system(cmd_process_4)


def Sub_figure(dir_path, args):
	isoform_proportion_inf = args.isoform_proportion_inf
	isoform_cpm_inf = args.isoform_cpm_inf
	group_info_inf = args.group_info_inf
	required_trans_inf = args.required_trans_inf
	bedgraph = args.bedgraph
	outf_dir = args.outf_dir.rstrip('/')
	figures = args.figures
	CDS_inf = args.CDS_inf
	genome_version = args.genome_version
	intron_shrinkage = args.intron_shrinkage
	order = args.order
	ref_gtf = args.ref_gtf
	input_gtf = args.input_gtf
	# 2.1 only show top 5 isoform's proportion and reshape the file format
	cmd_figure_1 = f"python {dir_path}/scripts/2_1_merge_isoforms_and_reshape_format.py {isoform_proportion_inf} {isoform_cpm_inf} {group_info_inf} {required_trans_inf} {outf_dir}"
	print(cmd_figure_1)
	os.system(cmd_figure_1)

	prop_reshaped_inf_name = f'{outf_dir}/' + re.sub(".txt", "_reshaped_merge_others.txt", isoform_proportion_inf.split('/')[-1])
	exp_reshaped_inf_name = f'{outf_dir}/' + re.sub(".esp|.txt", "_reshaped_merge_others.txt", isoform_cpm_inf.split('/')[-1])
	sorted_group_list = []
	with open(group_info_inf,'r') as inf:
		for index, line in enumerate(inf):
			if index == 0: continue
			arr = line.strip().split("\t")
			sorted_group_list.append(arr[0])
	sorted_group_info = ','.join(sorted_group_list)
	# 2.2 Generate bargraph and isoform structure figures [one example]
	if not os.path.exists(f"{outf_dir}/Example_res"):
		os.system(f"mkdir {outf_dir}/Example_res")
	if genome_version in ['hg19','GRCh37']:
		cmd_figure_2 = f"python {dir_path}/scripts/2_2_Generate_bar_structure_figure_example.py --gene [Ensembl_Gene_ID] --gene_name [Gene_Symbol] --transcript [Interested_Transcript_ID] --abundance_CPM_original {isoform_cpm_inf} --abundance_proportion {prop_reshaped_inf_name} --abundance_CPM {exp_reshaped_inf_name} --bedgraph {bedgraph} --sorted_group {sorted_group_info} --out_dir {outf_dir}/Example_res --figures {' '.join(figures)} --canonical_transcript {dir_path}/scripts/references/Gencode_v39_canonical_isoform.txt --basic_transcript {dir_path}/scripts/references/gencode.v34lift37.annotation_basic_trans.txt --anno_gtf {ref_gtf} --novel_gtf {input_gtf} --CDS_inf {CDS_inf} --genome_version {genome_version} --intron_shrinkage {intron_shrinkage} --order {order}"
	elif genome_version in ['hg38', 'GRCh38']:
		cmd_figure_2 = f"python {dir_path}/scripts/2_2_Generate_bar_structure_figure_example.py --gene [Ensembl_Gene_ID] --gene_name [Gene_Symbol] --transcript [Interested_Transcript_ID] --abundance_CPM_original {isoform_cpm_inf} --abundance_proportion {prop_reshaped_inf_name} --abundance_CPM {exp_reshaped_inf_name} --bedgraph {bedgraph} --sorted_group {sorted_group_info} --out_dir {outf_dir}/Example_res --figures {' '.join(figures)} --canonical_transcript {dir_path}/scripts/references/Gencode_v39_canonical_isoform.txt --basic_transcript {dir_path}/scripts/references/gencode.v39.annotation_basic_trans.txt --anno_gtf {ref_gtf} --novel_gtf {input_gtf} --CDS_inf {CDS_inf} --genome_version {genome_version} --intron_shrinkage {intron_shrinkage} --order {order}"
	with open(f"{outf_dir}/Template_to_generate_figures.sh", "w") as outf_figure_2:
		outf_figure_2.write(cmd_figure_2+'\n')


def Sub_difftest(dir_path, args):
	isoform_cpm_inf = args.isoform_cpm_inf
	tumor_num = args.tumor_num
	enriched_test_p = args.enriched_test_p
	enriched_test_tumor_cpm = args.enriched_test_tumor_cpm
	enriched_test_fc = args.enriched_test_fc
	#specificity_test_p = args.specificity_test_p
	specificity_test_tumor_cpm = args.specificity_test_tumor_cpm
	specificity_test_tissue_cpm = args.specificity_test_tissue_cpm
	specificity_test_tumor_percentage = args.specificity_test_tumor_percentage
	specificity_test_tissue_percentage = args.specificity_test_tissue_percentage
	specificity_test_exclude_tissue = args.specificity_test_exclude_tissue
	outf_dir = args.outf_dir.rstrip('/')
	## DE test ##
	cmd_difftest_1 = f"python {dir_path}/scripts/3_1_Tumor_vs_normal_DE_test.py --CPM_inf {isoform_cpm_inf} --Tumor_num {tumor_num} --cutoff_p {enriched_test_p} --cutoff_tumor_cpm {enriched_test_tumor_cpm} --cutoff_fc {enriched_test_fc} --outf_dir {outf_dir}"
	print(cmd_difftest_1)
	os.system(cmd_difftest_1)
	## prevalence test ##
	cmd_difftest_2 = f"python {dir_path}/scripts/3_2_Tumor_vs_normal_prevalence.py --CPM_inf {isoform_cpm_inf} --Tumor_num {tumor_num} --cutoff_tumor_cpm {specificity_test_tumor_cpm} --cutoff_tissue_cpm {specificity_test_tissue_cpm} --cutoff_tumor_percentage {specificity_test_tumor_percentage} --cutoff_tissue_percentage {specificity_test_tissue_percentage} --exclude_tissue {specificity_test_exclude_tissue} --outf_dir {outf_dir}"
	print(cmd_difftest_2)
	os.system(cmd_difftest_2)


def Sub_translation(dir_path, args):  # sourcery skip: extract-duplicate-method
	#genome = args.genome
	if args.genome_version in ['hg19','GRCh37','grch37']:
		genome = f'{dir_path}/scripts/references/hg19.fa'
	elif args.genome_version in ['hg38','GRCh38','grch38']:
		genome = f'{dir_path}/scripts/references/GRCh38.primary_assembly.genome.fa'
	input_trans_gtf = args.trans_gtf
	ref_gtf = args.ref_gtf
	outf_dir = args.outf_dir.rstrip('/')
	out_file = args.out_file
	isoform_cpm_inf = args.isoform_cpm_inf
	mode = args.mode

	#### 0. process reference gtf file to generate CDS annotation ####
	command_0 = (f"python {dir_path}/scripts/4_1_process_gtf.py {ref_gtf} {outf_dir}")
	print (command_0)
	os.system(command_0)
	print ('Step 0 completed: reference gtf has been converted to CDS annotation file.\n')

	#### 1. convert input gtf file into bed file ####
	command_1 = f"python {dir_path}/scripts/4_2_convert_gtf2transcript_1_to_1_based.py {input_trans_gtf} {outf_dir}/4_2_{out_file}_gtf_processed.txt"
	print (command_1)
	os.system(command_1)
	print ('Step 1 completed: input gtf has been converted to bed file.\n')

	#### 2. translate transcripts into proteins ####
	command_2 = f"python {dir_path}/scripts/4_3_seq_translate_coordinate.py --mode {mode} --trans_inf {outf_dir}/4_2_{out_file}_gtf_processed.txt --abundance_inf {isoform_cpm_inf} --gtf_inf {input_trans_gtf} --outf_name {outf_dir}/4_3_{out_file} --genome {genome} --cds {outf_dir}/4_1_Gencode_converted_CDS_annotation.txt"
	print (command_2)
	os.system(command_2)
	print ('Step 2 completed: input transcripts have been translated into proteins.\n')

	#### 3. classify transcripts and rename ####
	command_3 = f"python {dir_path}/scripts/4_4_classify_transcript.py {outf_dir}/4_3_{out_file}_protein.txt {outf_dir} {ref_gtf} {input_trans_gtf} {out_file}"
	print (command_3)
	os.system(command_3)
	print ('Step 3 completed: translated proteins have been classified into different categories.\n')


def Sub_CAR_T(dir_path, args):
	genome = args.genome_version
	tmhmm_dir = args.tmhmm_dir.rstrip('/')
	isoform_cpm_inf = args.isoform_cpm_inf
	protein_inf = args.protein_inf
	de_trans_inf = args.de_trans_inf
	spe_trans_inf = args.spe_trans_inf
	out_file = args.out_file
	tumor_num = args.tumor_num
	specificity_score = args.specificity_score
	tissue_cpm = args.tissue_cpm
	tissue_percentage = args.tissue_percentage
	isoform_proportion_inf = args.isoform_proportion_inf
	annotated_isoform_contri_inf = args.annotated_isoform_contri_inf
	trans_CDS_inf = args.trans_CDS_inf
	window_size = args.window_size
	ref_gtf = args.ref_gtf
	gencode_fasta = args.gencode_fasta	
	#other_isoform_contri_inf = args.other_isoform_contri_inf
	outf_dir = f"{args.outf_dir.rstrip('/')}/CAR_T"
	if not os.path.exists(outf_dir):
		os.system(f"mkdir {outf_dir}")

	#### 1. Perform TMHMM prediction based on protein sequence ####
	command_1 = f"{tmhmm_dir}/tmhmm {protein_inf} -short > {outf_dir}/5_1_{out_file}_tmhmm_analysis.txt"
	print (command_1)
	if not os.path.exists(f"{outf_dir}/5_1_{out_file}_tmhmm_analysis.txt"):
		os.system(command_1)
	print ('Step 1 completed: tmhmm prediction is done\n')

	#### 2. inference based on UniProt annotation ####
	command_2 = f"python {dir_path}/scripts/5_2_uniprot_annotation_infer.py {isoform_cpm_inf} {protein_inf} {genome} {gencode_fasta} {outf_dir} {out_file}"
	print (command_2)
	os.system(command_2)
	print ('Step 2 completed: inference based on UniProt annotation is done.\n')

	#### 3. combine results ####
	command_3 = (f"python {dir_path}/scripts/5_3_parse_result.py {out_file} {outf_dir}")
	print (command_3)
	os.system(command_3)
	print ('Step 3 completed: result is combined by considering both methods.\n')

	##### 4. integrate DE test and cell surface protein prediction ####
	command_4_a = f"python {dir_path}/scripts/5_4_Tumor_vs_normal_CAR_T_summary.py --outf_dir {outf_dir} --outf_key high_confidence --protein_inf {protein_inf} --cell_surface_inf {outf_dir}/5_3_{out_file}_high_confidence.txt --de_inf {de_trans_inf} --spe_inf {spe_trans_inf}"
	print (command_4_a)
	os.system(command_4_a)
	command_4_b = f"python {dir_path}/scripts/5_4_Tumor_vs_normal_CAR_T_summary.py --outf_dir {outf_dir} --outf_key less_confidence --protein_inf {protein_inf} --cell_surface_inf {outf_dir}/5_3_{out_file}_less_confidence.txt --de_inf {de_trans_inf} --spe_inf {spe_trans_inf}"
	print (command_4_b)
	os.system(command_4_b)
	print ('Step 4 completed: cell surface proteins that pass DE/prevalence test are obtained.\n')

	##### 5. prioritize CAR-T targets ####
	command_5_pre = f"python {dir_path}/scripts/5_5_prioritize_pre.py {tumor_num} {tissue_cpm} {tissue_percentage} {isoform_cpm_inf} {isoform_proportion_inf} {annotated_isoform_contri_inf} {outf_dir}/5_4_Tumor_vs_normal_CAR_T_high_confidence.txt {ref_gtf} {outf_dir}"
	print (command_5_pre)
	os.system(command_5_pre)

	command_6_post = f"python {dir_path}/scripts/5_6_prioritize_post.py {tumor_num} {specificity_score} {isoform_cpm_inf} {protein_inf} {trans_CDS_inf} {window_size} {outf_dir}"
	print (command_6_post)
	os.system(command_6_post)
	print ('Step 5 completed: prioritized CAR-T targets are obtained.\n')


def Sub_TCR(dir_path, args):
	genome = args.genome_version
	netMHCpan_dir = args.netMHCpan_dir.rstrip('/')
	HLA_str_inf = args.HLA_str_inf
	protein_inf = args.protein_inf
	isoform_cpm_inf = args.isoform_cpm_inf
	isoform_proportion_inf = args.isoform_proportion_inf
	annotated_isoform_contri_inf = args.annotated_isoform_contri_inf
	trans_CDS_inf = args.trans_CDS_inf
	de_trans_inf = args.de_trans_inf
	spe_trans_inf = args.spe_trans_inf
	tumor_num = args.tumor_num
	specificity_score = args.specificity_score
	binding_affi = args.binding_affi
	tissue_cpm = args.tissue_cpm
	tissue_percentage = args.tissue_percentage
	window_size = args.window_size
	out_file = args.out_file
	ref_gtf = args.ref_gtf
	outf_dir = f"{args.outf_dir.rstrip('/')}/TCR"
	if not os.path.exists(outf_dir):
		os.system(f"mkdir {outf_dir}")

	HLA_list = []
	with open(HLA_str_inf,"r") as HLA_inf:
		for index,line in enumerate(HLA_inf):
			if index == 0: continue
			arr = line.strip().split("\t")
			for each_arr in arr[1].split(","):
				if each_arr not in HLA_list:
					HLA_list.append(each_arr)
	HLA_str = ",".join(HLA_list)
	command_6_1 = f"python {dir_path}/scripts/6_1_perform_netMHCpan.py {netMHCpan_dir} {protein_inf} {HLA_str} {window_size} {outf_dir} {de_trans_inf} {spe_trans_inf}"
	print (command_6_1)
	if not os.path.exists(f"{outf_dir}/6_1_TCR_temp_out.txt"):
		os.system(command_6_1)
	print ('Step 1 completed: running netMHCpan.\n')

	command_6_2 = f"python {dir_path}/scripts/6_2_calculate_specificity_score_TCR.py {outf_dir}/6_1_Summarized_TCR_netMHCpan_reshaped_out.txt {protein_inf} {window_size} {HLA_str} {isoform_cpm_inf} {tumor_num} {trans_CDS_inf} {outf_dir}"
	print (command_6_2)
	os.system(command_6_2)
	print ('Step 2 completed: calculating tumor-specificity score for each peptides.\n')

	command_6_3 = f"python {dir_path}/scripts/6_3_final_prioritize.py {tumor_num} {specificity_score} {tissue_cpm} {tissue_percentage} {binding_affi} {isoform_proportion_inf} {annotated_isoform_contri_inf} {window_size} {ref_gtf} {outf_dir}"
	print (command_6_3)
	os.system(command_6_3)
	print ('Step 3 completed: prioritized TCR targets are obtained.\n')


def Sub_specificity(dir_path, args):
	transcript_ID = args.transcript_ID
	protein_inf = args.protein_inf
	tumor_num = args.tumor_num
	isoform_cpm_inf = args.isoform_cpm_inf
	cell_surface_inf = args.cell_surface_inf
	window_size = args.window_size
	start_site = args.start_site
	end_site = args.end_site
	outf_dir = f"{args.outf_dir.rstrip('/')}/Example_Specificity_{transcript_ID}"

	if not os.path.exists(outf_dir):
		os.system(f"mkdir {outf_dir}")

	#### 1. The matrix of peptide corresponding isoforms has been constructed ####
	command_1 = f"python {dir_path}/scripts/7_1_construct_peptide_matrix.py {transcript_ID} {protein_inf} {outf_dir} {window_size}"
	print (command_1)
	os.system(command_1)

	command_2 = f"python {dir_path}/scripts/7_2_construct_extracellular_matrix.py {cell_surface_inf} {protein_inf} {outf_dir}"
	print (command_2)
	os.system(command_2)
	print ('Step 1&2 completed: The matrix of peptide corresponding isoforms has been constructed.\n')

	#### 2. Calculation of tumor-specificity score ####
	command_3 = f"python {dir_path}/scripts/7_3_calculate_score_matrix.py {transcript_ID} {isoform_cpm_inf} {tumor_num} {outf_dir} {window_size}"
	print (command_3)
	os.system(command_3)
	print ('Step 3 completed: Calculation of tumor-specificity score.\n')

	#### 2.1 Read the position of the region with the highest specificity score #####
	region_inf_name = f"{outf_dir}/7_3_highest_score_region_WindowSize_{window_size}_{transcript_ID}.txt"
	total_window_size = 50
	with open(region_inf_name, 'r') as inf:
		for line in inf:
			arr = line.strip().split('\t')
			if arr[0] == 'Median':
				pos = int(float(arr[1].split(',')[0]))
				if (int(start_site) == 0) and (int(end_site) == 999999):
					start_site = max(1, pos - int(total_window_size/2))
					end_site = start_site + total_window_size - 1

	#### 3. The heatmap figure is generated ####
	command_3 = f"Rscript {dir_path}/scripts/7_4_heatmap_score_line.R {transcript_ID} {outf_dir} {window_size} {start_site} {end_site}"
	print (command_3)
	os.system(command_3)
	print ('Step 3 completed: The heatmap figure based on CPM is generated.\n')

	command_4 = f"Rscript {dir_path}/scripts/7_4_heatmap_score_line_z_score.R {transcript_ID} {outf_dir} {window_size} {start_site} {end_site}"
	print (command_4)
	os.system(command_4)
	print ('Step 4 completed: The heatmap figure based on z-score is generated.\n')

	#### 4. Exp line figure is generated ####
	command_5 = f"python {dir_path}/scripts/7_5_barplot_reshape.py {transcript_ID} {outf_dir} {window_size} {start_site} {end_site} {tumor_num}"
	print (command_5)
	os.system(command_5)

	command_6 = f"Rscript {dir_path}/scripts/7_6_barplot.R {transcript_ID} {outf_dir} {window_size} {pos}"
	print (command_6)
	os.system(command_6)
	print ('Step 5&6 completed: The box figure is generated.\n')


def Sub_topology(dir_path, args):
	transcript_ID = args.transcript_ID
	protein_inf = args.protein_inf
	gene_symbol = args.gene_symbol
	score_cutoff = float(args.score_cutoff)
	outf_dir = f"{args.outf_dir.rstrip('/')}/CAR_T/"
	command_1 = f"python {dir_path}/scripts/8_1_make_protter_figure.py {transcript_ID} {gene_symbol} {protein_inf} {outf_dir} {score_cutoff}"
	print (command_1)
	os.system(command_1)
	print ('Step 1 completed: The Protter figure is generated.\n')


if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write('[INFO] User interrupted; program terminated.')
		sys.exit(0)

