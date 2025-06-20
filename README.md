# `HEAL: Homeologous Exchange Automated Labeller`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

HEAL is a snakemake workflow for automated labelling of homeologous exchanges (HE). 

HEAL takes whole genome sequencing reads from allopolyploids and the genome assemblies and annotations of the two subgenomes/progenitors to infer the copy number over the genomes. It then uses synteny information to match the copy numbers of both subgenomes. The combined copy numbers provide valuable information regarding homeologous exchanges and other rearangements. 

## Workflow

The workflow can be summarized as follows:
- Reads are aligned to the appropriate subgenome and classified using [snake-EAGLE-RC](https://github.com/kenji-yt/snake-EAGLE-RC) (a snakemake workflow around [EAGLE-RC](https://github.com/tony-kuo/eagle?tab=readme-ov-file#eagle-rc)). For details see the README of snake-EAGLE-RC.
- Syntenic anchors (high confidence collinear homeologs), which define regions of synteny and therefore recombination candidates, are identified using [GENESPACE](https://github.com/jtlovell/GENESPACE) (by importing the module [snake-GENESPACE](https://github.com/kenji-yt/snake-GENESPACE) which is a snakemake workflow facilitating GENESPACE analysis).
- Bed files defining bins within the progenitor assemblies are created using [“bedtools makewindow”](https://github.com/arq5x/bedtools2). GC content is counted in these bins using “bedtools nuc”. Average per base mappability is computed in these bins using per base mappability obtained with [genmap](https://github.com/cpockrandt/genmap).
- Classified bam files, synteny information and binned mappability and GC statistics are then fed into the R-package [healr](https://github.com/kenji-yt/healr). This R-package provides allopolyploid species functions to infer and combine copy numbers between subgenomes. It also has functions for plotting and summary statistics. In HEAL it is used to:
    - Count reads in bins using ([featureCounts](https://subread.sourceforge.net/featureCounts.html)).
    - Filter and normalize based on mappability and GC content.
    - Infer copy number using the Circular Binary Segmentation algorithm implemented in [DNAcopy](https://rdrr.io/bioc/DNAcopy/).
    - Assign a copy number to each synthenic anchor based on the bin it overlaps with. 
    - Plot the combined copy number along each subgenomes gene order. 
    - Compute and return genome wide and per-chromosome statistics regarding the distribution of anchor sets (groups of corresponding anchors from different subgenomes) in each copy number ratio. 

## Installation 

- [install Snakemake via Conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
- git clone HEAL: 
```
git clone https://github.com/kenji-yt/HEAL.git
```

## Input 

HEAL requires an input directory formated in the following specific manner. First, the input directory should be named after the type of whole-genome sequencing (WGS) you have. The two options are "DNA" (for regular WGS) and "WGBS" (for whole-genome bisulfite sequencing). In this input directory you must have two directories called "progenitors" and "polyploids":
- The 'progenitors' directory should contain one subdirectory for each reference genome (each progenitor for allopolyploids). The name of each subdirectory should be the species or genome name and it will appear as such in the output files. Within each subdirectory you should have the corresponding reference genome in fasta/fastq format and annotation in gff format. Name these files as you wish as long as they have one of the following extensions: 
    - "fa","fasta","fq","fna","fastq" for the assemblies.
    - "gff", "gff3" for the annotations. 
- The 'polyploids' directory should contain one subdirectory per sample. The name of each subdirectory should be a unique sample identifier (don't give the same name to different sample directories). Each of these sample directories should only contain sequencing reads in fasta/fastq format (gzipped or not). There should be one file if the data is single-end and two files if paired-end. If paired-end, make sure the filenames are identical expect for a '1' and '2' indicating which side of the pair each file contains. All samples must have the same sequencing experiment type, so either all are paired-end or all are single-end. HEAL also assumes that the expected read length is the same for each sample. If you have samples with different read lengths and experiment type you must run HEAL separately for each.  


Your input directory should have the following structure:
```
input_directory/
├── progenitors/
│   ├── species_1
│   │   ├── annotation.gff
│   │   └── assembly.fa
│   └── species_2
|       ├── annotation.gff
│       └── assembly.fa
└── polyploids/
    ├── sample_1
    │   ├── reads_pe_1.fastq
    │   └── reads_pe_2.fastq
    ├── sample_2
    │   ├── reads_pe_1.fastq
    │   └── reads_pe_2.fastq
    └── sample_3
        ├── reads_pe_1.fastq
        └── reads_pe_2.fastq
```

HEAL can work with up to three subgenomes (hexaploid). The chokepoint is EAGLE-RC, which means if you can classify reads in another manner, then there is no more limit to the number of subgenomes. 

## Usage

You are now ready to run HEAL. From within the "HEAL/" directory downloaded from github, run the following command:
```
snakemake --use-conda --cores N --config INPUT_DIR='your/input/directory' 
```
Make sure to have snakemake make installed, to replace 'your/input/directory' with the path to your input directory and 'N' with the number of cores you wish to allocate to the job. If you installed snakemake in a conda environment make sure to activate it (eg. "conda activate snakemake_env").  

The outputs will now be generated in a results directory within the HEAL directory. 

#### Settings:

There are a number of settings you can modify for the HEAL analysis.
- **BIN_SIZE**: The size of bins in base pairs on which you infer copy number. The default value is 10000.
- **SOFT_CLIP**: If you wish to preserve reads that are softclipped in the bismark alignment, potentially to recover exact recombination breakpoints, you can do so by setting the config argument SOFT_CLIP to 'true' (set to 'False' by default). This will keep soft clipped reads in the alignment step. In our experience, these softclipped reads did not enable recovery of exact breakpoint. For DNA seq data, softclipped reads are always allowed by bwa-mem2. 
- **FILTER**: Whether of not to filter and trim the input reads using [fastp](https://github.com/OpenGene/fastp). You can set it to 'True', 'False' or to custom parameters to pass to fastp. Default value is 'True'. More details in the README of [snake-EAGLE-RC](https://github.com/kenji-yt/snake-EAGLE-RC). 
- **SAVE_HEALR_LISTS**: A boolean indicating if you wish to save the intermediate data of healr (healr formatted R-list as a directory). This is usefull if you want to manually explore the read counts (recommended) without having to run featureCounts again. Default is "True". 

So if you wish to run HEAL with 12 cores, 50kbp bins, soft-clipping allowed in bismark and no trimming and filtering:
```
snakemake --use-conda --cores 12 --config INPUT_DIR='/path/to/WGBS' SOFT_CLIP='True' BIN_SIZE=50000 FILTER='False' 
```
#### Bin Size

If you realize that the bin size is not appropriate and want to try another value, you do not need to repeat all the GENESPACE, alignment and read classification steps. To repeat the analysis with a different bin size, move the reproducibility report and the healr directory out of the results. Then, run the snakemake command with the new BIN_SIZE value. 

## Output 

Results will be written to a directory called "results" inside the HEAL directory. In this directory you will find the following files and directories:
    - eagle_rc: Contains the eagle installation and one directory per sample with the eagle-rc results and a script used to produce these results.
    - fastp (if FILTER='True'): Contains one directory per sample with filtered and trimmed read files and quality check reports. 
    - qualimap: Contains one directory per sample containing the output of qualimap for every ".ref." bam files in the 'eagle_rc' directory. 
    - bismark/bwa: Contains bam files of the reads aligned to each subgenome. 
    - genespace: Contains a directory for MCScanX, a directory for scripts and 'run_dir' which contains the formated input data for GENESPACE as well as the results of the GENESPACE analysis. For more details see [GENESPACE](https://github.com/jtlovell/GENESPACE).
    - genmap: Contains the output of genmap for each subgenome. 
    - healr: Contains an input directory for healr ('input_dir') which can be used for further healr analyses. It also contains one directory called stats which has data tables  and one dire
    - logs: Contains logs for each analysis.
    - MultiQC: Contains the file "multiqc_report.html" which compiles qualimap and fastp reports.  
    - snakemake_EAGLE_RC_reproducibility_report.txt: A text file with details about the input and output files and the tools and parameters used. 
