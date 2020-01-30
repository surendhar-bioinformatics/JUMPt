#Contents of this file
Introduction
Release notes
Software Requirements
Hardware Requirements
Installation
Command Line Arguments
Maintainers

# Introduction
JUMPt (JUMP-turnover) software tool is for detrmining the protein turnover rates in metabolically labeled animals using mass spectrometry(MS) data. JUMPt uses novel differential equation-based mathematical model to determine the reliable and accurate protein turnover rates. The proposed method determines the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously. JUMPt is part of JUMP Software Suite (shortly JUMP), which is an integrative omics data processing and analysis tool, including protein/peptide database creation, database search, identification filtering, quantification, network analysis, proteogenomics and protein turnover analysis .

#Release notes
Version 1:

In this version 
1. We assume that the overall amount of proteins in mice is unchanged over time as the mice used in this study are adult. 
2. The Lys used in protein synthesis originate from the food intake, with the rest recycled through protein degradation. 
3. oncentration of soluble Lys is conserved, the Lys that is absorbed from the food excreted by degradation pathways. 

#Software Requirements
The program is written by a combination of Python3 and MATLAB. It should run on every system with Python3 with required modules loaded and MATLAB R2017 and above. The minimum required Perl version should be Perl 5.8 or better.

Perl modules are needed:

Parallel::ForkManager
Statistics::Distributions
Clone
Excel/Writer/XLSX.pm (optional but recommended)
To install perl modules, please refer to: http://www.cpan.org/modules/INSTALL.html


#Hardware Requirements
Starting from JUMPg_v2.3.1, the program can be run on either high performance computing systme or a single server.

The program has been successfully tested on the following system:

Cluster mode (key parameters: 'cluster = 1' & 'Job_Management_System = SGE'):

Batch-queuing system: SGE, version 6.1u5, qsub and qstat
128 GB memory and 64 cores on each node
Single server mode (key parameters: 'cluster = 0' & 'processors_used = 8'):

32 GB memory
2 GHz CPU processors with 8 cores
Installation
After download the source code, you can put it in any working directory (e.g. /home/usr/JUMPg). IMPORTANT: The folder containing all the source code needs to be accessible by each node.

#INSTALLATION:

step 0: unzip the source code and test data packages in the current directory by running the following commands:

unzip JUMPg_v2.3.1.zip (already unzipped if downloaded from github)
unzip exampleData_v6.2.zip
Step 1: set up module and program paths by running the following commands:

cd programs
perl path_setup.pl
cd ..
Step 2: download AnnoVar annotation files by running the following commands (take human hg19 as an example):

mkdir annotations
cd annotations
perl ../programs/c/customizedDB/annovar/annotate_variation.pl -downdb -buildver hg19 -downdb knownGene human/
cd ..
Step 3: prepare reference genome (FASTA format) # this file is required for splice junction and RNA 6FT analysis

For model organisms, the reference genome should be available through UCSC genome browser database (http://genome.ucsc.edu/).

Example command: wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz Concatenate all the fasta files using the following command cat *.fa > hg19_genome_raw.fa

Step 4: (Optional but recommended) update uniProt-UCSC ID convertion table based on pairwise sequence alignment.

For human, an updated uniProt-UCSC ID convertion table (named 'ID_convert_uniProt2UCSChg19.txt') is available in the exampleData_v6.2/ folder. This table is used to correct the UCSC downloaded table, using the following command:
perl tools/updateIdTab.pl annotations/human/hg19_kgXref.txt exampleData_v6.2/resources/ID_convert_uniProt2UCSChg19.txt
Input File Preparation
(You can skip this section if you only want to run through the example.)

Besides MS data, the program takes up to three types of genomics/NGS data as inputs for the customized peptide database building.

For building the mutation peptide database: an AnnoVar input format file that contains the genomic locations and allele (ref/alt) information file is required.
For building the junction peptide database: a splice junction file reported by STAR (or with equivalent format) is required.
For building the RNA 6FT database: RNA nucleotide sequences in FASTQ format is required (can be either de novo assemebled transcripts or reference mRNA sequences).
Users can find example input files in exampleData_v6.2.tar.gz (within the rna_database/ folder).

If users process the data from RNA-seq raw reads, to assist the input file preparation process, a wrapper (called 'RNAseq_preprocess.pl' located in the tools/ folder) is provided that contains RNAseq alignment (by STAR), mutation detection (by GATK), and novel splice junction filtering.

To run the wrapper, the following 3rd party software is required:

STAR: https://github.com/alexdobin/STAR/releases Note that the reference genome should be index by the STAR 'genomeGenerate' module
Picard: http://broadinstitute.github.io/picard/
GATK: https://www.broadinstitute.org/gatk/ Note that GATK requires the reference genome to be sorted by kayrotypic order
Samtools: http://samtools.sourceforge.net/ Only used for indexing BAM files
Users can check each requirement by filling the parameter file of the wrapper ('RNAseq_preprocess.params' located in the tools/ folder).

Once finished with the steps above, the wrapper can be executed using the following command:

perl tools/RNAseq_preprocess_v1.0.0.pl RNAseq_preprocess.params <test_PE1.fq> <test_PE2.fq> <test_PE1.fq> <test_PE2.fq> Pair 1 and 2 of RNAseq FASTAQ files

The wrapper will output:

Mutation file in AnnoVar input format: filtered_output.annovar
Novel splice junction file in STAR format: SJ.out.filtered2uniq.tab
Users can take these two files as input for the JUMPg program.

Run the example
Since the program supports multistage database search analysis, we have designed a two-stage test dataset. For the 1st stage, the MS/MS spectra are searched against a peptide database pooled from uniProt, mutations and splice junctions; the matching results are filtered to ~1% false discovery rate. For the 2nd stage, the remaining high quality spectra are searched against the 6FT database.

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

#Maintainers
To submit bug reports and feature suggestions, please contact:

surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

#Acknowledgement
We gratefully acknowledge St. Jude Children’s Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite. 
