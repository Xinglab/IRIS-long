import snakemake.utils

snakemake.utils.min_version('6.5.3')

configfile: 'snakemake_config.yaml'

onsuccess:
	print('workflow success')

onerror:
	print('workflow error')

DEFAULT_MEM_MB=20 * 1024  # 20 GB
DEFAULT_TIME_HOURS=24

# Specifying this as an input to a rule will disable that rule.
# This can be used in combination with "ruleorder:" to determine what
# rule should be used to create a particular output file.
UNSATISFIABLE_INPUT='unsatisfiable_input_file_path'

def all_input(wildcards):
	inputs = dict()
	inputs.update(check_all_results())
	return inputs

localrules: all
rule all:
	input:
		unpack(all_input),

def check_all_results():
	outputs = dict()
	outputs['CAR_T_result'] = os.path.join(config['outf_dir'], "CAR_T", "5_5_Summarized_CAR_T_prioritized_targets_final.txt")
	outputs['TCR_result'] = os.path.join(config['outf_dir'], "TCR", "6_3_Summarized_TCR_prioritized_targets.txt")
	outputs['Figure_sh'] = os.path.join(config['outf_dir'], "Template_to_generate_figures.sh")
	return outputs

def get_genome_fasta(wildcards):
	if (config['genome_fasta_path']) and (config['genome_fasta_path'] != ''):
		return config['genome_fasta_path']
	else:
		if config['genome_version'] in ['hg38','GRCh38']:
			return os.path.join(config['IRIS_long_path'], 'scripts', 'references', 'GRCh38.primary_assembly.genome.fa')
		elif config['genome_version'] in ['hg19','GRCh37']:
			return os.path.join(config['IRIS_long_path'], 'scripts', 'references', 'hg19.fa')

def get_genome_gtf(wildcards):
	if (config['genome_gtf_path']) and (config['genome_gtf_path'] != ''):
		return config['genome_gtf_path']
	else:
		if config['genome_version'] in ['hg38','GRCh38']:
			return os.path.join(config['IRIS_long_path'], 'scripts', 'references', f"genome.v{config['gencode_release_version']}.annotation.gtf")
		elif config['genome_version'] in ['hg19','GRCh37']:
			return os.path.join(config['IRIS_long_path'], 'scripts', 'references', 'genome.v34lift37.annotation.gtf')

def get_gencode_fasta(wildcards):
	if config['genome_version'] in ['hg38','GRCh38']:
		return os.path.join(config['IRIS_long_path'], 'scripts', 'references', f"gencode.v{config['gencode_release_version']}.pc_translations.fa")
	elif config['genome_version'] in ['hg19','GRCh37']:
		return os.path.join(config['IRIS_long_path'], 'scripts', 'references', 'gencode.v34lift37.pc_translations.fa')

rule download_genome_fasta:
	output:
		genome_fasta=get_genome_fasta(None),
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'download_genome_fasta.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'download_genome_fasta.err'),
	params:
		genome_fasta_url=config['genome_fasta_url'],
	shell:
		'curl -L {params.genome_fasta_url} -o {output.genome_fasta}.gz && '
		'gunzip {output.genome_fasta}.gz 1>> {log.out} 2>> {log.err}'

rule download_genome_gtf:
	output:
		genome_gtf=get_genome_gtf(None),
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'download_genome_gtf.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'download_genome_gtf.err'),
	params:
		genome_gtf_url=config['genome_gtf_url'],
	shell:
		'curl -L {params.genome_gtf_url} -o {output.genome_gtf}.gz && '
		'gunzip {output.genome_gtf}.gz 1>> {log.out} 2>> {log.err}'

rule download_additional_file:
	output:
		gencode_fasta=get_gencode_fasta(None),
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'download_additional_file.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'download_additional_file.err'),
	params:
		gencode_fasta_url=f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{config['gencode_release_version']}/gencode.v{config['gencode_release_version']}.pc_translations.fa.gz",
	shell:
		'curl -L {params.gencode_fasta_url} -o {output.gencode_fasta}.gz && '
		'gunzip {output.gencode_fasta}.gz 1>> {log.out} 2>> {log.err}'

rule pre_processing:
	input:
		transcript_gtf = config['input_transcript_gtf'],
		transcript_abundance = config['input_transcript_abundance_matrix'],
		download_gencode_fasta = get_gencode_fasta(None),
		download_genome_gtf = get_genome_gtf(None),
		download_genome_fasta = get_genome_fasta(None),
	output:
		trans_CPM = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM.txt"),
		trans_proportion = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_proportion.txt"),
		anno_trans_contribution = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_gene_annotated_isoform_contribution.txt"),
		bed_file = os.path.join(config['outf_dir'], "samples_BedGraph.bed"),
	params:
		conda_wrapper = config['conda_wrapper'],
		IRIS_long = os.path.join(config['IRIS_long_path'], 'IRIS_long_main.py'),
		normalized_mode = config['normalized_mode'],
		ref_gtf = get_genome_gtf,
		outf_dir = config['outf_dir'],
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'Pre_processing.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'Pre_processing.err'),
	resources:
		mem_mb=config['general_mem_gb'] * 1024,
		time_hours=config['general_time_hr'],
	shell:
		'{params.conda_wrapper} python {params.IRIS_long} Preprocess' 
		' --input_gtf {input.transcript_gtf}'
		' --input_abundance {input.transcript_abundance}'
		' --normalized_mode {params.normalized_mode}'
		' --ref_gtf {params.ref_gtf}'
		' --outf_dir {params.outf_dir}'
		' 1> {log.out}'
		' 2> {log.err}'

rule DE_transcripts:
	input:
		trans_CPM = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM.txt"),
	output:
		TE_trans = os.path.join(config['outf_dir'], "3_1_Tumor_vs_normal_DE_test.txt"),
		TS_trans = os.path.join(config['outf_dir'], "3_2_Tumor_vs_normal_prevalence.txt"),
	params:
		conda_wrapper = config['conda_wrapper'],
		IRIS_long = os.path.join(config['IRIS_long_path'], 'IRIS_long_main.py'),
		outf_dir = config['outf_dir'],
		num_tumor_sample = config['num_tumor_sample'], 
		enriched_test_p_value = config['enriched_test_p_value'], 
		enriched_test_cpm_in_tumor = config['enriched_test_cpm_in_tumor'], 
		enriched_test_fc = config['enriched_test_fc'],
		specificity_test_cpm_in_tumor = config['specificity_test_cpm_in_tumor'],
		specificity_test_cpm_in_tissue = config['specificity_test_cpm_in_tissue'],
		specificity_test_minimum_percentage_in_tumor = config['specificity_test_minimum_percentage_in_tumor'],
		specificity_test_maximal_percentage_in_tissue = config['specificity_test_maximal_percentage_in_tissue'],
		specificity_test_exclude_tissue = config['specificity_test_exclude_tissue'],
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'DE_transcripts.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'DE_transcripts.err'),
	resources:
		mem_mb=config['general_mem_gb'] * 1024,
		time_hours=config['general_time_hr'],
	shell:
		'{params.conda_wrapper} python {params.IRIS_long} DiffTest' 
		' --isoform_cpm_inf {input.trans_CPM}'
		' --tumor_num {params.num_tumor_sample}'
		' --enriched_test_p {params.enriched_test_p_value}'
		' --enriched_test_tumor_cpm {params.enriched_test_cpm_in_tumor}'
		' --enriched_test_fc {params.enriched_test_fc}'
		' --specificity_test_tumor_cpm {params.specificity_test_cpm_in_tumor}'
		' --specificity_test_tissue_cpm {params.specificity_test_cpm_in_tissue}'
		' --specificity_test_tumor_percentage {params.specificity_test_minimum_percentage_in_tumor}'
		' --specificity_test_tissue_percentage {params.specificity_test_maximal_percentage_in_tissue}'
		' --specificity_test_exclude_tissue {params.specificity_test_exclude_tissue}'
		' --outf_dir {params.outf_dir}'
		' 1> {log.out}'
		' 2> {log.err}'

rule translation:
	input:
		transcript_gtf = config['input_transcript_gtf'],
		trans_CPM = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM.txt"),
	output:
		PC_fasta = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_PC.fasta"),
		PC_NMD_fasta = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_PC_and_NMD.fasta"),
		CDS_info_table = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_detailed_match_ID.txt"),
	params:
		conda_wrapper = config['conda_wrapper'],
		IRIS_long = os.path.join(config['IRIS_long_path'], 'IRIS_long_main.py'),
		outf_dir = config['outf_dir'],
		mode = config['mode'],
		ref_gtf = get_genome_gtf,
		tumor_name = config['tumor_key_word'],
		genome_version = config['genome_version'],
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'Translation.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'Translation.err'),
	resources:
		mem_mb=config['general_mem_gb'] * 1024,
		time_hours=config['general_time_hr'],
	shell:
		'{params.conda_wrapper} python {params.IRIS_long} Translation' 
		' --mode {params.mode}'
		' --trans_gtf {input.transcript_gtf}'
		' --isoform_cpm_inf {input.trans_CPM}'
		' --genome_version {params.genome_version}'
		' --ref_gtf {params.ref_gtf}'
		' --out_file {params.tumor_name}'
		' --outf_dir {params.outf_dir}'
		' 1> {log.out}'
		' 2> {log.err}'

rule CAR_T_prediction:
	input:
		PC_fasta = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_PC.fasta"),
		trans_CPM = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM.txt"),
		trans_proportion = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_proportion.txt"),
		anno_trans_contribution = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_gene_annotated_isoform_contribution.txt"),
		CDS_info_table = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_detailed_match_ID.txt"),
		TE_trans = os.path.join(config['outf_dir'], "3_1_Tumor_vs_normal_DE_test.txt"),
		TS_trans = os.path.join(config['outf_dir'], "3_2_Tumor_vs_normal_prevalence.txt"),
	output:
		CAR_T_output= os.path.join(config['outf_dir'], "CAR_T", "5_5_Summarized_CAR_T_prioritized_targets_final.txt"),
	params:
		conda_wrapper = config['conda_wrapper'],
		IRIS_long = os.path.join(config['IRIS_long_path'], 'IRIS_long_main.py'),
		outf_dir = config['outf_dir'],
		tmhmm_dir = config['tmhmm_dir'],
		num_tumor_sample = config['num_tumor_sample'],
		tumor_name = config['tumor_key_word'],
		ref_gtf = get_genome_gtf,
		ref_gencode_fasta = get_gencode_fasta,
		genome_version = config['genome_version'],
		minimum_specificity_score_CAR_T = config['minimum_specificity_score_CAR_T'],
		maximal_peptide_cpm_in_tissue_CAR_T = config['maximal_peptide_cpm_in_tissue_CAR_T'],
		maximal_tissue_percentage_CAR_T = config['maximal_tissue_percentage_CAR_T'],
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'CAR_T.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'CAR_T.err'),
	resources:
		mem_mb=config['general_mem_gb'] * 1024,
		time_hours=config['general_time_hr'],
	shell:
		'{params.conda_wrapper} python {params.IRIS_long} CAR_T' 
		' --tmhmm_dir {params.tmhmm_dir}'
		' --tumor_num {params.num_tumor_sample}'
		' --protein_inf {input.PC_fasta}'
		' --isoform_cpm_inf {input.trans_CPM}'
		' --isoform_proportion_inf {input.trans_proportion}'
		' --de_trans_inf {input.TE_trans}'
		' --spe_trans_inf {input.TS_trans}'
		' --genome_version {params.genome_version}'
		' --annotated_isoform_contri_inf {input.anno_trans_contribution}'
		' --out_file {params.tumor_name}'
		' --specificity_score {params.minimum_specificity_score_CAR_T}'
		' --tissue_cpm {params.maximal_peptide_cpm_in_tissue_CAR_T}'
		' --tissue_percentage {params.maximal_tissue_percentage_CAR_T}'
		' --trans_CDS_inf {input.CDS_info_table}'
		' --ref_gtf {params.ref_gtf}'
		' --gencode_fasta {params.ref_gencode_fasta}'
		' --outf_dir {params.outf_dir}'
		' 1> {log.out}'
		' 2> {log.err}'

rule TCR_prediction:
	input:
		PC_NMD_fasta = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_PC_and_NMD.fasta"),
		trans_CPM = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM.txt"),
		trans_proportion = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_proportion.txt"),
		anno_trans_contribution = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_gene_annotated_isoform_contribution.txt"),
		CDS_info_table = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_detailed_match_ID.txt"),
		TE_trans = os.path.join(config['outf_dir'], "3_1_Tumor_vs_normal_DE_test.txt"),
		TS_trans = os.path.join(config['outf_dir'], "3_2_Tumor_vs_normal_prevalence.txt"),
	output:
		TCR_output= os.path.join(config['outf_dir'], "TCR", "6_3_Summarized_TCR_prioritized_targets.txt"),
	params:
		conda_wrapper = config['conda_wrapper'],
		IRIS_long = os.path.join(config['IRIS_long_path'], 'IRIS_long_main.py'),
		outf_dir = config['outf_dir'],
		netMHCpan_dir = config['netMHCpan_dir'],
		num_tumor_sample = config['num_tumor_sample'],
		tumor_name = config['tumor_key_word'],
		ref_gtf = get_genome_gtf,
		genome_version = config['genome_version'],
		minimum_specificity_score_TCR = config['minimum_specificity_score_TCR'],
		maximal_peptide_cpm_in_tissue_TCR = config['maximal_peptide_cpm_in_tissue_TCR'],
		maximal_tissue_percentage_TCR = config['maximal_tissue_percentage_TCR'],
		HLA_str_inf = config['HLA_str_inf'],
		AA_window_size = config['AA_window_size'],
		binding_affi = config['binding_affi'],
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'TCR.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'TCR.err'),
	resources:
		mem_mb=config['general_mem_gb'] * 1024,
		time_hours=config['general_time_hr'],
	shell:
		'{params.conda_wrapper} python {params.IRIS_long} TCR' 
		' --netMHCpan_dir {params.netMHCpan_dir}'
		' --HLA_str_inf {params.HLA_str_inf}'
		' --isoform_cpm_inf {input.trans_CPM}'
		' --tumor_num {params.num_tumor_sample}'
		' --protein_inf {input.PC_NMD_fasta}'
		' --de_trans_inf {input.TE_trans}'
		' --spe_trans_inf {input.TS_trans}'
		' --window_size {params.AA_window_size}'
		' --isoform_proportion_inf {input.trans_proportion}'
		' --genome_version {params.genome_version}'
		' --annotated_isoform_contri_inf {input.anno_trans_contribution}'
		' --out_file {params.tumor_name}'
		' --binding_affi {params.binding_affi}'
		' --specificity_score {params.minimum_specificity_score_TCR}'
		' --tissue_cpm {params.maximal_peptide_cpm_in_tissue_TCR}'
		' --tissue_percentage {params.maximal_tissue_percentage_TCR}'
		' --trans_CDS_inf {input.CDS_info_table}'
		' --ref_gtf {params.ref_gtf}'
		' --outf_dir {params.outf_dir}'
		' 1> {log.out}'
		' 2> {log.err}'

rule figure:
	input:
		trans_CPM = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM.txt"),
		trans_proportion = os.path.join(config['outf_dir'], "samples_abundance_combined_CPM_proportion.txt"),
		group_info_inf = config["group_info_inf"],
		required_trans_inf = config["required_trans_inf"],
		bed_file = os.path.join(config['outf_dir'], "samples_BedGraph.bed"),
		CDS_info_table = os.path.join(config['outf_dir'], f"4_4_{config['tumor_key_word']}_detailed_match_ID.txt"),
	output:
		Figure_template = os.path.join(config['outf_dir'], "Template_to_generate_figures.sh"),
	params:
		conda_wrapper = config['conda_wrapper'],
		IRIS_long = os.path.join(config['IRIS_long_path'], 'IRIS_long_main.py'),
		ref_gtf = get_genome_gtf,
		transcript_gtf = config['input_transcript_gtf'],
		outf_dir = config['outf_dir'],
		genome_version = config['genome_version'],
	log:
		out=os.path.join(config['outf_dir'], 'log_dir', 'Figure.out'),
		err=os.path.join(config['outf_dir'], 'log_dir', 'Figure.err'),
	resources:
		mem_mb=config['general_mem_gb'] * 1024,
		time_hours=config['general_time_hr'],
	shell:
		'{params.conda_wrapper} python {params.IRIS_long} Figure' 
		' --isoform_proportion_inf {input.trans_proportion}'
		' --isoform_cpm_inf {input.trans_CPM}'
		' --group_info_inf {input.group_info_inf}'
		' --required_trans_inf {input.required_trans_inf}'
		' --bedgraph {input.bed_file}'
		' --outf_dir {params.outf_dir}'
		' --CDS_inf {input.CDS_info_table}'
		' --figures Isoform Single_isoform Structure'
		' --genome_version {params.genome_version}'
		' --ref_gtf {params.ref_gtf}'
		' --input_gtf {params.transcript_gtf}'
		' 1> {log.out}'
		' 2> {log.err}'
