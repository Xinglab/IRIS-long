# IRIS-long: Isoform peptides from RNA splicing for Immunotherapy target Screening - based on Long-read RNA-seq data

## Table of Contents

* [Overview](#overview)
* [Dependencies](#dependencies)
* [Usage](#usage)
  + [Data processing](#data-processing)
  + [Data visualization and analysis](#data-visualization-and-analysis)


## Overview

For the quantification part using ESPRESSO, please refer to [ESPRESSO GitHub page](https://github.com/Xinglab/espresso) and [TEQUILA-seq GitHub page](https://github.com/Xinglab/TEQUILA-seq).

<img src="./files/TEQUILA-seq_Analysis_Workflow.png" width="800"/>

Our data processing scripts are designed to work with raw Oxford Nanopore (ONT) signal data (FAST5 format) as input. However, these scripts can work with data from any long-read sequencing platform (e.g., PacBio) as long as the data is provided in FASTQ/SAM/BAM format. These scripts encompass the following steps:
1. **Basecalling**: Raw ONT signal data (FAST5 format) is basecalled into nucleotide sequences (FASTQ format) using Guppy in fast mode.
2. **Alignment**: Basecalled reads are mapped to a user-supplied reference genome using [minimap2](https://github.com/lh3/minimap2) (Li H., *Bioinformatics* 2018) together with user-supplied reference transcript annotations.
3. **Transcript isoform discovery and quantification**: Full-length transcript isoforms are discovered and quantified from long-read RNA-seq alignment files using [ESPRESSO](https://github.com/Xinglab/espresso).

We also have scripts designed to visualize and further characterize transcript isoforms identified from TEQUILA-seq data. These scripts can perform the following tasks:
1. **Visualize discovered transcript isoforms**: Given a collection of samples subjected to TEQUILA-seq, we can visualize the structures and relative abundances of all transcript isoforms discovered for a given gene.
2. **Detect group-specific and sample-specific transcript isoforms**: For a collection of samples subjected to TEQUILA-seq, if the samples can be partitioned into different groups, we can identify transcript isoforms with group-specific expression and usage. Similarly, we can also identify transcript isoforms with sample-specific expression and usage.
3. **Characterize alternative splicing events underlying discovered transcript isoforms**: Local differences in transcript structure between a given isoform and the canonical isoform of the corresponding gene are classified into different alternative splicing patterns, including exon skipping, alternative 5' and 3' splice sites, mutually exclusive exons, retained introns, alternative first or last exons, and complex splicing. 
4. **Predict NMD-targeted transcript isoforms**: Transcript isoforms targeted by mRNA nonsense-mediated decay (NMD) are predicted from the set of isoforms discovered from TEQUILA-seq data based on the 50 nt rule.

## Dependencies

To run our scripts, the following dependencies will need to be installed and available on `$PATH`:

* [Snakemake](https://snakemake.readthedocs.io) (v5.31.1) 
* Guppy (must be downloaded manually from the [ONT software download page](https://community.nanoporetech.com/downloads) since a login is required
* [minimap2](https://github.com/lh3/minimap2) (v2.17)
* [SAMtools](http://samtools.sourceforge.net) (v1.9)
* [BLAST](https://www.ncbi.nlm.nih.gov/blast/) (v2.10.1)
* [HMMER](http://hmmer.org/) (v3.3.1)
* [UCSC KentUtils](http://hgdownload.soe.ucsc.edu/admin/exe/)
  + bedGraphToBigWig
  + faToTwoBit
  + twoBitInfo
* [Python](https://www.python.org/) 3.8
  + [NumPy](https://numpy.org/) (v1.20.1)
  + [pandas](https://pandas.pydata.org/) (v1.1.4)
  + [SciPy](https://scipy.org/) (v1.5.4)
  + [statsmodels](https://www.statsmodels.org/) (v0.12.2)
  + [NetworkX](https://networkx.org/) (v2.6.3)
  + [BeautifulSoup4](https://pypi.org/project/beautifulsoup4/) (v4.8.2)
  + [ConfigArgParse](https://pypi.org/project/ConfigArgParse/)
* [R](https://www.r-project.org/) (v4.0.5)
  + [ggplot2](https://ggplot2.tidyverse.org/)
  + [tidyverse](https://www.tidyverse.org/)
  + [ggplotify](https://cran.r-project.org/package=ggplotify)
  + [scales](https://scales.r-lib.org/)
  + [forcats](https://forcats.tidyverse.org/)
  + Check for thread support with `perl -e 'use threads; print("ok\n")'`

## Usage

### Combine ESPRESSO result (optional)

A sub-command is used to combine ESPRESSO results from different runs. The command can start from raw ESPRESSO gtf and ESPRESSO abundance matrix and it will generate combined gtf file and combined isoform abundance matrix.

Our script can be run as follows:

```
python /mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/IRIS_long/IRIS_long_main.py Combine [-h] --gtf_list /path/to/espresso_gtf_file/list --outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --gtf_list                                          Path to espresso_gtf_file list

    --outf_dir                                          Folder of output 

```

`gtf_list` lists ESPRESSO gtf files from different runs. 

An example `gtf_list` file would be:
```
./samples_N2_R0_updated.gtf    PARCB
/mnt/isilon/xing_lab/aspera/xuy/CloneTechTissueAll_ESPRESSO_0225/samples_N2_R0_updated_hg38.gtf    Tissue
```
Note: columns are separated by `tab`. 


### Data processing

A sub-command is used to pre-process ESPRESSO results. The command can start from raw ESPRESSO gtf and ESPRESSO abundance matrix and it will generate bed format file (derived from gtf file) as well as normalized abundance matrix using CPM as unit in both transcript level and gene level. Besides, it will also generate isoform proportion matrix, which would be used in the following steps

Our script can be run as follows:

```
python /mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/IRIS_long/IRIS_long_main.py Preprocess [-h] --espresso_gtf /path/to/espresso_gtf_file --espresso_abundance /path/to/espresso_abundance_matrix_file --normalized_mode /choose/from/'SAM'/or/'ESPRESSO' --folder_sam /path/to/folder/of/sam_files --outf_dir /path/to/folder/of/output/file

script arguments:
    -h, --help                                          Show this message and exit

    --espresso_gtf                                      ESPRESSO gtf file

    --espresso_abundance                                ESPRESSO abundance file

    --normalized_mode                                   Choose to normalize CPM based on SAM files or ESPRESSO output file

    --folder_sam                                        Folder of sam files

    --outf_dir                                          Folder of output 

```

### Data visualization and analysis

A sub-command is used to visualize process results. The command can start from results generated from previous step, and it will finally generate bar-graph figures for both isoform proportion and isoform abundance (CPM) in a gene, as well as the transcript structure figures for all involved isoforms in a gene.

Our script can be run as follows:

```
python /mnt/isilon/xing_lab/aspera/xuy/snakemake_ESPRESSO_reference/pipeline_test/IRIS_long/IRIS_long_main.py Figure --isoform_porportion_inf /path/to/isoform/proportion/matrix --isoform_cpm_inf /path/to/isoform/abundance/matrix/CPM --group_info_inf /path/to/file/containing/group_info --required_trans_inf /path/to/file/containing/required_transcripts --bedgraph /path/to/processed/bed/file --outf_dir /path/to/folder/of/output/file --figures Isoform Single_isoform Structure

script arguments:
    -h, --help                                          Show this message and exit

    --isoform_porportion_inf                            Isoform proportion matrix from previous step

    --isoform_cpm_inf                                   Isoform abundance matrix file (CPM) from previous step

    --group_info_inf                                    Sample group information

    --required_trans_inf                                Transcripts need to show in final figure

    --bedgraph                                          Bed file generated from previous step

    --outf_dir                                          Folder of output

    --figures                                           Figures expected to generate, could be multiple choices from ['Isoform','Single_isoform','Structure'], seperated by ' ' (white space)

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
Thus, we need to input all the interested transcripts in  `required_trans_inf` so that they could be included in the genrated figures. 

An example `required_trans_inf` file would be:
```
Gene_ID Trans_ID  Gene_symbol
ENSG00000026508 ENST00000434472;ENST00000428726 CD44
```
Note: columns are separated by `tab`. Multiple required transcripts are separated by `;`.

