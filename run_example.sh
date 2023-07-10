## IRIS-long on a toy dataset (chromsome X on melanoma + tissue samples)

## Data processing
python ./IRIS_long_main.py Preprocess --espresso_gtf example_dataset/samples_updated_combined.gtf --espresso_abundance example_dataset/samples_abundance_combined.txt --folder_sam . --normalized_mode ESPRESSO --outf_dir ./example_dataset/

## Differential test
python ./IRIS_long_main.py DiffTest --isoform_cpm_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO.txt --tumor_num 49 --detest_p 0.05 --detest_tumor_cpm 3 --detest_fc 2 --pretest_p 1e-6 --pretest_tumor_cpm 3 --pretest_tissue_cpm 1 --outf_dir ./example_dataset/

## Transcript Translation
python ./IRIS_long_main.py Translation --mode long-read --trans_gtf example_dataset/samples_updated_combined.gtf --isoform_cpm_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO.txt --genome_version hg38 --ref_gtf scripts/references/gencode.v39.annotation.gtf --out_file Melanoma --outf_dir ./example_dataset/

## CAR-T target prediction
python ./IRIS_long_main.py CAR_T --tmhmm_dir /mnt/isilon/xing_lab/aspera/xuy/tmhmm-2.0c/bin/ --tumor_num 49 --protein_inf example_dataset/4_4_Melanoma_PC.fasta --isoform_cpm_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO.txt --isoform_proportion_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO_proportion.txt --genome_version hg38 --annotated_isoform_contri_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO_gene_annotated_isoform_contribution.txt --out_file Melanoma --specificity_score 2 --outf_dir ./example_dataset/ --tissue_cpm 20 --tissue_number 3 --trans_CDS_inf example_dataset/4_4_Melanoma_detailed_match_ID.txt

## TCR target prediction
#python ./IRIS_long_main.py TCR  --netMHCpan_dir /home/xuy2/xuy2/program/netMHCpan-4.1 --HLA_str HLA-A02:01 --isoform_cpm_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO.txt --tumor_num 49 --protein_inf example_dataset/4_4_Melanoma_PC.fasta --window_size 9 --isoform_proportion_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO_proportion.txt --genome_version hg38 --annotated_isoform_contri_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO_gene_annotated_isoform_contribution.txt --outf_dir ./example_dataset/ --specificity_score 3 --binding_affi 500 --tissue_cpm 10 --tissue_number 3 --trans_CDS_inf example_dataset/4_4_Melanoma_detailed_match_ID.txt

## Example visualization
python ./IRIS_long_main.py Figure --isoform_proportion_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO_proportion.txt --isoform_cpm_inf example_dataset/samples_abundance_combined_CPM_ESPRESSO.txt --group_info_inf example_dataset/group_info.txt --required_trans_inf example_dataset/required_trans.txt --bedgraph example_dataset/samples_BedGraph.bed --outf_dir ./example_dataset/ --figures Isoform Single_isoform Structure --genome_version hg38
# Input infomration of target transcript in ./example_dataset/Template_to_generate_figures.sh
bash ./example_dataset/Template_to_generate_figures_L1CAM.sh
# Generated figures could be found in ./example_dataset/Example_res

