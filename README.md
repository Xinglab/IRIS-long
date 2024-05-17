## IRIS-long: Isoform peptides from RNA splicing for Immunotherapy target Screening - based on Long-read RNA-seq data

## Table of Contents

* [Overview](#overview)
* [Installation](#installation)
* [Dependencies](#dependencies)
* [Usage (Snakemake)](#usage-in-snakemake-recommended)
* [Usage (Command)](#usage-in-command)
* [Example](#example)

## Overview

IRIS-long tool is designed to discover novel tumor antigen from RNA dysregulation for immunotherapy, and it works with transcript-level quantification result based on long-read RNA-seq data. If you start with fast5 raw files, please refer to [ESPRESSO GitHub page](https://github.com/Xinglab/espresso) and [TEQUILA-seq GitHub page](https://github.com/Xinglab/TEQUILA-seq) for the transcript identification and quantification; in which, Guppy (Basecalling), [minimap2](https://github.com/lh3/minimap2) (Alignment) and [ESPRESSO](https://github.com/Xinglab/espresso) (Quantification) tool might be used.

<img src="./files/IRIS_long_workflow_diagram.png" width="800"/>

## Installation

The IRIS-long program can be downloaded directly from the repository, and when `install` commond is run, a conda environment will be created. Detailed commands are shown below:

```
git clone https://github.com/Xinglab/IRIS-long.git
cd ~/IRIS-long
chmod +x install
./install
```

IRIS-long requires annotated genome fasta and gtf files, user-specified files should be downloaded firstly and moved into `~/IRIS_long/IRIS_long/scripts/references/`, otherwise genome fasta (hg38) and corresponding genocode gtf (genome.v39.annotation.gtf) will be downloaded and used automatically. 


## Dependencies

To run our scripts, the following dependencies will need to be installed and available on `$PATH`:

* [SAMtools](http://samtools.sourceforge.net) 
* [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)
* [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) 
* [Python3](https://www.python.org/)
  + [NumPy](https://numpy.org/) 
  + [ConfigArgParse](https://pypi.org/project/ConfigArgParse/) 
  + [SciPy](https://scipy.org/) 
  + [statsmodels](https://www.statsmodels.org/) 
  + [BioPython](https://biopython.org/) 
  + [Matplotlib](https://matplotlib.org/)
* [R](https://www.r-project.org/) 
  + [ggplot2](https://ggplot2.tidyverse.org/)
  + [tidyverse](https://www.tidyverse.org/)
  + [ggplotify](https://cran.r-project.org/package=ggplotify)
  + [cowplot](https://github.com/wilkelab/cowplot)
  + [scales](https://scales.r-lib.org/)
  + [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
  + [viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)


## Usage in Snakemake (Recommended)

The most straightforward method to execute the IRIS-long tool is by utilizing the built-in Snakemake workflow. Please modify the corresponding parameters in the [snakemake_config.yaml](./snakemake_config.yaml) file, as well as [group_info.txt](./group_info.txt) and [required_trans.txt](./required_trans.txt) files for figure generation step.

To run IRIS-long, simply submit this command:
```
sbatch ./run
```

To generate the directed acyclic graph (DAG) plot of whole IRIS-long workflow, we could run
```
./conda_wrapper snakemake --profile ./snakemake_profile --dag | dot -Tsvg > IRIS_long_dag.svg
```

For detailed explanation of other parameters, please refer to corresponding section in [Usage in Command](#usage-in-command).



## Usage in Command

**Important Note**: 
1. The input parameter for `--outf_dir` in all steps below should be the same folder, which should be the folder contains ESPRESSO output abundance matrix and gtf file.
2. The columns in all files should be separated by single `tab` rather than `white space`.

### Data integration (Optional, not included in Snakemake)

This sub-command is used to combine ESPRESSO results from different runs. The command can start from raw ESPRESSO gtf and ESPRESSO abundance matrix and it will generate combined gtf file and combined isoform abundance matrix.

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py Combine [-h] \
--allowed_dist /allowed/distance/for/each/ends/to/collapse/novel/transcripts \
--gtf_list /path/to/espresso_gtf_file/list \
--outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --gtf_list                                          Path to espresso_gtf_file list

    --outf_dir                                          Folder of output 

```

`gtf_list` lists ESPRESSO gtf files from different runs. 
Each row of `gtf_list` should contain `<gtf_file> <abundance_matrix> <group_name>`

An example `gtf_list` file would be:
```
./tumor.gtf    ./tumor_abundance.esp    Tumor
./tissue_gtf   ./tissue_abundance.esp   Tissue
```
Note: columns are separated by `tab`. 



### Data processing

**It is important to note that when running the IRIS-long tool, the order of samples should be sorted by tumor and normal tissues. Specifically, the columns corresponding to tumor samples should be placed ahead (on the left) of the columns representing normal tissues in the expression matrix file.**

This sub-command is used to pre-process ESPRESSO results. The command can start from raw ESPRESSO gtf and ESPRESSO abundance matrix and it will generate bed format file (derived from gtf file) as well as normalized abundance matrix using CPM as unit in both transcript level and gene level. Besides, it will also generate isoform proportion matrix, which would be used in the following steps

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py Preprocess [-h] \
--espresso_gtf /path/to/espresso_gtf_file \
--espresso_abundance /path/to/espresso_abundance_matrix_file \
--normalized_mode /choose/from/'SAM'/or/'ESPRESSO' \
--ref_gtf /gencode/gtf/file \
--folder_sam /path/to/folder/of/sam_files \
--outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --espresso_gtf                                      ESPRESSO gtf file

    --espresso_abundance                                ESPRESSO abundance file

    --ref_gtf                                           Reference Gencode gtf

    --normalized_mode                                   Choose to normalize CPM based on SAM files or ESPRESSO output file, choose from ['SAM','ESPRESSO']

    --folder_sam                                        Folder of sam files

    --outf_dir                                          Folder of output 

```



### Differential test

This sub-command is used to perform differential tests between tumor samples and normal tissue samples. This step consists of two test (based on isoform level): differential expression test (two-sided wilcoxon-ranksum test) and prevalence test (fisher-exact test).

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py DiffTest [-h] \
--isoform_cpm_inf /path/to/isoform_cpm_matrix \
--tumor_num /number/of/tumor/samples \
--enriched_test_p /p-value/in/tumor-enriched test \
--enriched_test_tumor_cpm /Cutoff/of/CPM/across/tumor/in/tumor-enriched test \
--enriched_test_fc /Cutoff/of/fold/change/in/tumor-enriched test \
--specificity_test_tumor_cpm /Cutoff/of/CPM/across/tumor/in/tumor-specificity test \
--specificity_test_tumor_percentage /Minimum/precentage/of/tumor/samples/that/express/given/transcript \
--specificity_test_tissue_cpm /Cutoff/of/CPM/across/tissue/in/tumor-specificity test \
--specificity_test_tissue_percentage /Maximum/precentage/of/tissue/samples/that/express/given/transcript \
--specificity_test_exclude_tissue /Excluded/tissues/in/tumor-specificity/transcript/identification \
--outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --isoform_cpm_inf                                   Isoform CPM file, e.g. samples_abundance_combined_CPM_ESPRESSO.txt

    --tumor_num                                         Number of tumor samples

    --enriched_test_p                                   Cutoff of p-value in tumor-enriched test, default = 0.05

    --enriched_test_tumor_cpm                           Cutoff of median CPM value in tumor samples, used to decide whether an isoform is highly expressed, default = 5

    --enriched_test_fc                                  Cutoff of fold change between tumor and normal, used to decide DE isoform, default = 2

    --specificity_test_tumor_cpm                        Cutoff of CPM value in tumor samples, used to decide whether an isoform is considered as expressed (default = 5)

    --specificity_test_tumor_percentage                 Minimum precentage of tumor samples that express given transcript (default = 50%, which is 0.5)

    --specificity_test_tissue_cpm                       Cutoff of CPM value in tissue samples, used to decide whether an isoform is considered as expressed (default = 1)

    --specificity_test_tissue_percentage                Maximum precentage of tissue samples that express given transcript (default = 10%, which is 0.1)

    --specificity_test_exclude_tissue                   Excluded tissues in tumor-specificity transcript identification, separated by comma. (default = "Testis")

    --outf_dir                                          Folder of output 

```



### Transcript Translation

This sub-command is used to classify transcripts into different types (protein-coding, NMD or fragment). Then translate protein-coding transcripts into proteins.

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py Translation [-h] \
--mode /short-read/or/long-read \
--trans_gtf /path/to/ESPRESSO/gtf/file \
--isoform_cpm_inf /path/to/isoform_cpm_matrix \
--genome_version /hg19/or/hg38 \
--ref_gtf /path/to/reference/gencode/gtf/file \
--out_file /prefix/of/name/of/output/file \
--outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --mode                                              Long-read or short-read RNA-seq data mode, default is long-read

    --trans_gtf                                         Generated gtf file, e.g. samples_updated_combined.gtf

    --isoform_cpm_inf                                   Isoform CPM file, e.g. samples_abundance_combined_CPM_ESPRESSO.txt

    --genome_version                                    Choose from ['GRCH38','GRCH37','hg38','hg19']

    --ref_gtf                                           Reference gencode annotation, e.g. ./references/gencode.v39.annotation.gtf

    --out_file                                          Prefix of the name of output file

    --outf_dir                                          Folder of output 

```



### CAR-T target prediction

This sub-command is used to decide protein topology and further discover potential targets for CAR-T therapy.

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py CAR_T [-h] \
--tmhmm_dir /path/of/tmhmm \
--protein_inf /path/to/generated/protein/fasta \
--isoform_cpm_inf /path/to/isoform_cpm_matrix \
--isoform_proportion_inf /path/to/isoform/proportion/matrix \
--annotated_isoform_contri_inf /path/to/file \
--trans_CDS_inf /path/to/file \
--genome_version /hg19/or/hg38 \
--tumor_num /number/of/tumor/samples \
--specificity_score /cutoff/of/specificity_score \
--tissue_cpm /cutoff/of/transcripts/in/tissue/samples \
--tissue_percentage /maximum/tolerable/percentage/of/tissues/that/are/allowed/to/have/transcript/higher/than/the/given/CPM/expression/threshold \
--ref_gtf /gencode/gtf/file \
--gencode_fasta /reference/gencode/translation/fasta \
--de_trans_inf /Tumor-enriched/transcripts \
--spe_trans_inf /Tumor-specific/transcripts \
--window_size /Window/size \
--out_file /prefix/of/name/of/output/file \
--outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --tmhmm_dir                                         File path of TMHMM tool (directory is needed)

    --protein_inf                                       Generated protein fasta file

    --isoform_cpm_inf                                   Isoform CPM file, e.g. samples_abundance_combined_CPM_ESPRESSO.txt

    --isoform_proportion_inf                            Isoform proportion file, e.g. samples_abundance_combined_CPM_ESPRESSO_proportion.txt

    --annotated_isoform_contri_inf                      File generated before, which ends with "_annotated_isoform_contribution.txt"

    --trans_CDS_inf                                     File generated before, format of which is like "4_4_*_detailed_match_ID.txt"

    --genome_version                                    Choose from ['GRCH38','GRCH37','hg38','hg19']

    --tumor_num                                         Number of tumor samples

    --specificity_score                                 Cutoff of specificity score (default = 1)

    --tissue_cpm                                        Cutoff of (maximum tolerable) CPM of transcripts encode given peptide in tissue samples (default = 10)

    --tissue_percentage                                 Maximum tolerable percentage of tissues that are allowed to have transcript higher than the given CPM expression threshold (default = 20%, which is 0.2)

    --ref_gtf                                           Reference gencode annotation, e.g. gencode.v39.annotation.gtf

    --gencode_fasta                                     Reference gencode translation.fasta, e.g. gencode.v39.pc_translations.fa

    --de_trans_inf                                      Tumor-enriched transcripts, e.g. 3_1_Tumor_vs_normal_DE_test.txt

    --spe_trans_inf                                     Tumor-specific transcripts, e.g. 3_2_Tumor_vs_normal_prevalence.txt

    --window_size                                       Window size (default = 9 AAs)

    --out_file                                          Prefix of the name of output file

    --outf_dir                                          Folder of output 

```



### TCR target prediction

This sub-command is used to predict samples-specific HLA types and further discover potential targets for TCR therapy.

When peforming TCR target prediction job, we need to input HLA alleles as the parameter, which could be obtained based on bam/sam files by tools such as [HLA-LA](https://github.com/DiltheyLab/HLA-LA) for given samples, or we can manually specify comman HLA alleles such as HLA-A02:01,HLA-A01:01

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py TCR [-h] \
--netMHCpan_dir /path/of/netMHCpan \
--HLA_str_inf /HLA/alleles/inf \
--protein_inf /path/to/generated/protein/fasta \
--isoform_cpm_inf /path/to/isoform_cpm_matrix \
--isoform_proportion_inf /path/to/isoform/proportion/matrix \
--genome_version /hg19/or/hg38 \
--annotated_isoform_contri_inf /path/to/file \
--trans_CDS_inf /path/to/file \
--tumor_num /number/of/tumor/samples \
--specificity_score /cutoff/of/specificity_score \
--tissue_cpm /cutoff/of/transcripts/in/tissue/samples \
--tissue_percentage /maximum/tolerable/percentage/of/tissues/that/are/allowed/to/have peptide/higher/than/the/given/CPM/expression/threshold \
--binding_affi /cutoff/of/binding/affinity/between/HLA/complex/and/peptide \
--window_size /size/of/sliding/window \
--ref_gtf /gencode/gtf/file \
--de_trans_inf /Tumor-enriched/transcripts \
--spe_trans_inf /Tumor-specific/transcripts \
--out_file /prefix/of/name/of/output/file \
--outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --netMHCpan_dir                                     File path of netMHCpan tool (directory is needed)

    --HLA_str_inf                                       File containing HLA alleles information, first column is sample, and second column is interesed HLA allele that separated by comma

    --protein_inf                                       Generated protein fasta file, such as 4_4_XXX_PC.fasta

    --isoform_cpm_inf                                   Isoform CPM file, e.g. samples_abundance_combined_CPM_ESPRESSO.txt

    --isoform_proportion_inf                            Isoform proportion file, e.g. samples_abundance_combined_CPM_ESPRESSO_proportion.txt

    --genome_version                                    Choose from ['GRCH38','GRCH37','hg38','hg19']

    --annotated_isoform_contri_inf                      File generated before, which is like "_annotated_isoform_contribution.txt"

    --trans_CDS_inf                                     File generated before, which is like  "4_4_*_detailed_match_ID.txt"

    --tumor_num                                         Number of tumor samples

    --binding_affi                                      Cutoff of binding affinity between HLA complex and peptides (default = 500)

    --specificity_score                                 Cutoff of specificity score (default = 3)

    --tissue_cpm                                        Cutoff of (maximum tolerable) CPM of transcripts encode given peptide in tissue samples (default = 10)

    --tissue_percentage                                 Maximum tolerable percentage of tissues that are allowed to have transcript higher than the given CPM expression threshold (default = 20%, which is 0.2)

    --ref_gtf                                           Reference gencode annotation, e.g. gencode.v39.annotation.gtf

    --de_trans_inf                                      Tumor-enriched transcripts, e.g. 3_1_Tumor_vs_normal_DE_test.txt

    --spe_trans_inf                                     Tumor-specific transcripts, e.g. 3_2_Tumor_vs_normal_prevalence.txt

    --window_size                                       Window size (default = 9 AAs)

    --out_file                                          Prefix of the name of output file

    --outf_dir                                          Folder of output 

```



### Example visualization

The command based on the results generated from previous step, and it will generate a bash file `Template_to_generate_figures.sh` as the output. The bar-graph figures for both isoform proportion and isoform abundance (CPM) in a gene, as well as the transcript structure figure would be generated when interested `Ensembl_Gene_ID`, `Gene_Symbol` and `Ensembl_Transcript_ID` are specified.

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py Figure [-h] \
--isoform_proportion_inf /path/to/isoform/proportion/matrix \
--isoform_cpm_inf /path/to/isoform/abundance/matrix/CPM \
--group_info_inf /path/to/file/containing/group_info \
--required_trans_inf /path/to/file/containing/required_transcripts \
--bedgraph /path/to/processed/bed/file \
--genome_version /hg19/or/hg38 \
--CDS_inf /generated/CDS/file \
--intron_shrinkage /Intron/shrinkage/fold/in/isoform/structure/figure \
--ref_gtf /gencode/gtf/file \
--espresso_gtf /novel/gtf/file \
--outf_dir /path/to/folder/of/output/file \
--figures Isoform Single_isoform Structure

script arguments:
    -h, --help                                          Show this message and exit

    --isoform_proportion_inf                            Isoform proportion infile, e.g. samples_abundance_combined_CPM_ESPRESSO_proportion.txt

    --isoform_cpm_inf                                   Isoform CPM file, e.g. samples_abundance_combined_CPM_ESPRESSO.txt

    --group_info_inf                                    Sample group information

    --required_trans_inf                                Transcripts need to show in final figure

    --bedgraph                                          Generated bedgraph file for sample, e.g. samples_BedGraph.bed

    --genome_version                                    Choose from ['GRCH38','GRCH37','hg38','hg19']

    --CDS_inf                                           Generated CDS file, e.g. 4_4_*_detailed_match_ID.txt

    --intron_shrinkage                                  Intron shrinkage fold in isoform structure figure

    --ref_gtf                                           Reference Gencode gtf

    --espresso_gtf                                      ESPRESSO gtf file

    --outf_dir                                          Folder of output

    --figures                                           Figures expected to generate, could be multiple choices from ['Isoform','Single_isoform','Structure'], seperated by ' ' (white space, no quotation mark)

```

`group_info_inf` indicates how many sample groups we have, and how many samples in each group. 

An example `group_info_inf` file would be:
```
Group   Number_of_samples
Tumor   16
Tissue  30
```
Note: columns are separated by `tab`. And the order is important, based on this example, we know the first 16 samples in transcript expression matrix belong to Tumor group and the rest 30 samples are normal tissues.


`required_trans_inf` indicates what transcript we want show in the figure. In default, the five transcripts we would include are: 
1. The interested transcript
2. The canonical transcript in a gene (if it's not the interested transcript, otherwise it would be the longest annotated transcript, based on Gencode annotation)
3. The 3rd - 5th transcripts would be the transcripts with the highest average proportion across all samples among the rest transcripts in a gene. 
Thus, we need to input all the interested transcripts in  `required_trans_inf` so that they could be included in the genrated figures (one gene per row).

An example `required_trans_inf` file would be:
```
Gene_ID Trans_ID  Gene_symbol
ENSG00000026508 ENST00000434472;ENST00000428726 CD44
```
Note: columns are separated by `tab`. Multiple required transcripts are separated by `;`.



### Tumor specificity

The step is to calculate the tumor-specificity score for each region (e.g. 9 AAs) along the given transcript-derived protein sequence (from predicted CAR-T targets). Besides, it will also generate the figure showing the change of tumor-specificity scores along the protein sequence.

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py Specificity [-h] \
--transcript_ID /EnsemblID/of/interested/transcript \
--tumor_num /number/of/tumor/samples \
--protein_inf /path/to/generated/protein/fasta \
--isoform_cpm_inf /path/to/isoform/abundance/matrix/CPM \
--cell_surface_inf /predicted/cell/surface/protein/in/given/samples \
--window_size /size/of/sliding/window \
--start_site /starting/position/of/visualized/region \
--end_site /ending/position/of/visualized/region \
--outf_dir /path/to/folder/of/output/file 


script arguments:
    -h, --help                                          Show this message and exit

    --transcript_ID                                     Ensembl ID of interested transcript

    --tumor_num                                         Number of tumor samples

    --protein_inf                                       Generated protein fasta file, which is like "4_4_*_PC.fasta"

    --isoform_cpm_inf                                   Isoform CPM file, e.g. samples_abundance_combined_CPM_ESPRESSO.txt

    --cell_surface_inf                                  Generated cell surface proteins file, which is like "5_3_*_high_confidence.txt"

    --window_size                                       Size of sliding window, (default = 9)

    --start_site                                        Starting position of visualized protein region (default shows the 50 AAs region with the highest tumor-specificity)

    --end_site                                          Ending position of visualized protein region (default shows the 50 AAs region with the highest tumor-specificity)

    --outf_dir                                          Folder of output

```



### Protein topology 

The step is to generate protein topology figure using [Protter](https://wlab.ethz.ch/protter/help/) tool, and this step is only for predicted CAR-T targets, please run it after CAR-T prediction step.
Note, this step is using API service from Protter tool, there is no need to download the tool.

Our script can be run as follows:

```
python ~/IRIS_long/IRIS_long_main.py Topology [-h] \
--transcript_ID /EnsemblID/of/interested/transcript \
--gene_symbol /Interested/gene/name \
--protein_inf /path/to/generated/protein/fasta \
--score_cutoff /tumor/specificity/score/cutoff \
--outf_dir /path/to/folder/of/output/file 


script arguments:
    -h, --help                                          Show this message and exit

    --transcript_ID                                     Ensembl ID of interested transcript

    --gene_symbol                                       Gene name of interested gene

    --protein_inf                                       Generated protein fasta file, which is like "4_4_*_PC.fasta"

    --score_cutoff                                      Specificity score cutoff (default = 3)

    --outf_dir                                          Folder of output

```


### Additionally useful scripts (Only limited to Xing lab)

If interested target is due to differential gene expression between tumor and normal tissue group, we could try to query the expression profile of given gene from IRIS db in CHOP HPC. This script will generate a file containing TPM value of given gene across TCGA and GTEx samples, based on short-read RNA-seq data. Besides, it will also generate a box plot accordingly.

```
python ~/IRIS_long/scripts/Supp_2_extract_gene_expression.py [gene_symbol] [Output_dir] (Tumor_type) (Tissue_type)

Such as:
### Include all TCGA cancers and all GTEx normal tissues
python ~/IRIS_long/scripts/Supp_2_extract_gene_expression.py HIPK2  .       

### Only include TCGA-BRCA cancer and GTEx tissues that in ClonTech tissue list
python ~/IRIS_long/scripts/Supp_2_extract_gene_expression.py HIPK2  . BRCA ClonTech         
```


If interested splice junction is involved in classical alternative splicing events, we could try to query this splicing junction from IRIS db in CHOP HPC. This script will generate a file containing normalized read counts mapped to given splicing junction across TCGA and GTEx samples, based on short-read RNA-seq data. Besides, it will also generate a box plot accordingly.

```
python ~/IRIS_long/scripts/Supp_4_extract_SJC.py [gene_symbol] [SJ_coordinate (hg19)] [Output_dir]

Such as:

python ~/IRIS_long/scripts/Supp_4_extract_SJC.py HIPK2 chr7:139299240:139305146 .
```


If we want to validate specific junction/exon of interested isoform, we could extrac the reads that only mapped to given gene from given sample, then IGV tool could be used to the following visualization.

```
python ~/IRIS_long/scripts/Supp_3_extract_sam.py [gene_symbol] [sample] [genome_version] [Sam_folder] [Output_dir]

python ~/IRIS_long/scripts/Supp_3_select_reads_map_interested_junction.py [transcript_ID] [gene_symbol] [region_left_boundary] [region_right_boundary] [ESPRESSO_folder] [Output_dir]

Such as:

python ~/IRIS_long/scripts/Supp_3_extract_sam.py L1CAM M202 hg38 /home/xuy2/xuy2/Stored_scratch/snakemake_1.3.1_Melanoma/alignment .

python ~/IRIS_long/scripts/Supp_3_select_reads_map_interested_junction.py ENST00000370055 L1CAM 153872697 153875761 /home/xuy2/xuy2/Stored_scratch/snakemake_1.3.1_Melanoma/espresso_out/q_work_dir .
```


## Example
A toy dataset (gtf and expression matrix) has been inclued in `./example_dataset/`, which is about the identified transcripts from chromosome X across melanoma cell lines and normal tissues.

Please follow the instruction on `run_example.sh` under `IRIS_long` folder:
```
bash run_example.sh
```

Predicted potential CAR-T targets could be found in `./example_dataset/CAR_T/5_5_Summarized_CAR_T_prioritized_targets_final.txt`

Predicted potential TCR targets could be found in `./example_dataset/TCR/6_3_Summarized_TCR_prioritized_targets.txt`

After modifying [./example_dataset/Template_to_generate_figures.sh](./example_dataset/Template_to_generate_figures.sh) accordingly, the generated figures of input transcript could be found in [./example_dataset/Example_res/](./example_dataset/Example_res/). 

For example:
<img src="./example_dataset/Example_res/Figure_L1CAM_bed_list_sorted_ENST00000370055_rank_single_color.png" width="800"/>
<img src="./example_dataset/Example_res/Bar_sample_isoform_L1CAM_ENST00000370055_exp.png" width="800"/>

