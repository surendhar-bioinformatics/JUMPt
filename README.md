# Contents of this file
- Introduction
- Release notes
- Software Requirements
- Hardware Requirements
- Installation
- Command Line Arguments
- Maintainers
- Acknoledgements

# Introduction
JUMPt (JUMP-turnover) software tool is for detrmining the protein turnover rates in metabolically labeled animals using mass spectrometry(MS) data. JUMPt uses novel differential equation-based mathematical model to determine the reliable and accurate protein turnover rates. The proposed method determines the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously. 

JUMPt is part of JUMP Software Suite (shortly JUMP), which is an integrative omics data processing and analysis tool, including protein/peptide database creation, database search, identification filtering, quantification, network analysis, proteogenomics and protein turnover analysis .

# Release notes
Version 1:

In this version 
1. We assume that the overall amount of proteins in mice is unchanged over time as the mice used in this study are adult. 
2. The Lys used in protein synthesis originate from the food intake, with the rest recycled through protein degradation. 
3. oncentration of soluble Lys is conserved, the Lys that is absorbed from the food excreted by degradation pathways. 

# Software Requirements
The program is written in Python3 and MATLAB language. It should run on every system with Python3 with required modules loaded and MATLAB R2017 and above. The minimum required Perl version should be Perl 5.8 or better.

Python modules are needed:

The following python modules are required. For new user, we encourage to use Anacoda, a free and open-source distribution of the Python.
- numpy
- scipy
- pandas
- matplotlib
- subprocess
To install Anaconda distribution, please refer to: https://www.anaconda.com/

MATLAB toolbox needed:
We recommed to use MATLAB R2014 (The MathWorks, Inc., Natick, Massachusetts, United States) or above version on Linux, Mac or Windows. Our program is mainly use the following toolbox other than basic toolboxes.
- Global Optimization tool box


# Hardware Requirements
The program can be run on either linux or windows. Memory size of 8GB is an optimum and higher memory may required depending on size of the data.

The current program has been successfully tested on the following system:
16 GB memory
3.3 GHz CPU processors with 6 cores

# Installation
Installation of the script is not required. Download all the scripts to any working directory (e.g. /home/usr/JUMPt). IMPORTANT: All the scripts including associated modules (associted with the program) should be placed in the same folder. Once the scripts are saved, open and run "main.py" in python console.

# Run the example

A testing dataset (test_data.xlsx) is available along with the scripts for evaluation purposes. Start the "main.py" script in Python console. and follow the instruction. Input your data file (preferably in .xslx) along with full path. 
During the excution of program, expected output is saved in "Results" folder which is created in the wordking directory. Run time for demo on a desktop computer is about 2-10 min depending on the system configuretation.

The program perform multistep optimization to find the optimal degredation rates (turnover rates or Half-lives). The output of each step is saved in For the 1st stage, the MS/MS spectra are searched against a peptide database pooled from uniProt, mutations and splice junctions; the matching results are filtered to ~1% false discovery rate. For the 2nd stage, the remaining high quality spectra are searched against the 6FT database.

1st stage analysis:

Step 1: cp 'jump_g_v2.2.stage1.params' (included in exampleData_v6.2/parameters folder) to your working directory and edit the following parameters:

The following assumes that your working directory is at "/home/usr/JUMPg"

input_ref_proteins = /home/usr/JUMPg/exampleData_v6.2/uniProt/reference_uniProt_human.fasta
input_mutation = /home/usr/JUMPg/exampleData_v6.2/rna_database/mutations.txt
input_junction = /home/usr/JUMPg/exampleData_v6.2/rna_database/junctions.txt
annovar_annotation_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_knownGene.txt
gene_ID_convert_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_kgXref.txt
reference_genome = /home/usr/genomes/hg19_genome_raw.fa
Step 2: cp spectra file 'exampleData_v6.2/spectra/MS_test.mzXML' (included in exampleData_v6.2.tar.gz) to your working diredirectory and run the program using the following command:

perl /home/usr/JUMPg_v2.3.1/programs/JUMPg_v2.3.1.pl jump_g_v2.2.stage1.params MS_test.mzXML

Output: the program will generate a folder named 'gnm_round1_test1'. The final results are all included in the 'publications' folder that includes 6 files:

identified_peptide_table.txt: a text table, showing the identified peptide sequences, PSM counts, tags and scores of the best PSM, and database entries.
mutation_peptides.txt: a text table showing the identified mutation peptides
mutation_peptides.bed: mutation peptides with genomic locations in BED format, which can be co-displayed with other genomic information on the UCSC genome browser
junction_peptides.txt: a text table showing the identified novel junction peptides
junction_peptides.bed: novel junction peptides with genomic locations in BED format for visualization on the UCSC genome browser
reference_peptides.bed: reference peptides with genomic locations in BED format for visualization on the UCSC genome browser
For multistage analysis, the program also generates the unmatched high quality MS/MS spectra, of which the path is recorded in this file: gnm_round1_test1/multistage/qc_MSMS_input.txt. This file will be taken as input for 2nd stage analysis.

2nd stage analysis:

Step 1: cp 'jump_g_v2.2.stage2.params' (included in exampleData_v6.2/parameters) to your working directory and edit the following parameters:

input_transcript_6FT = /home/usr/JUMPg/exampleData_v6.2/rna_database/assembled_transcripts.fasta
annovar_annotation_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_knownGene.txt
gene_ID_convert_file = /home/usr/JUMPg/JUMPg_v2.3.1/annotations/human/hg19_kgXref.txt
reference_genome = /home/usr/genomes/hg19_genome_raw.fa
reference_protein = /home/usr/JUMPg/exampleData_v6.2/uniProt/reference_uniProt_human.fasta
Step 2: copy qc_MSMS_input.txt (that records the path to unmatched high quality MS/MS spectra) to current directory:

cp gnm_stage1_test1/multistage/qc_MSMS_input.txt .

Step 3: run the program by the command:

perl /home/usr/JUMPg_v2.3.1/programs/JUMPg_v2.3.1.pl jump_g_v2.2.stage2.params qc_MSMS_input.txt

Output: similar to the 1st stage result, the program will generate a folder named 'gnm_round2_test1' containing results in its 'publications' folder.

# Maintainers
To submit bug reports and feature suggestions, please contact:

surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

# Acknowledgement
We gratefully acknowledge St. Jude Children’s Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite. 
