# Contents of this file
- Introduction
- Release notes
- Software Requirements
- Hardware Requirements
- Installation
- Input File Preparation
- Run JUMPt program (test data provided) 
- Maintainers
- Acknowledgments

# Introduction
JUMPt (JUMP-turnover) software determines the protein turnover rates in metabolically labeled animals using mass spectrometry (MS) data. JUMPt uses a novel differential equation-based mathematical model to calculate the reliable and accurate protein turnover rates. The proposed method calculates the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously.

JUMPt is part of JUMP Software Suite (shortly JUMP), is an ongoing large software suite developed for the need of mass spectrometry (MS)- based proteomics, metabolomics, and the integration with genomics for network analysis at the level of systems biology. Currently, JUMP can handle protein/peptide database creation, database search, identification filtering, quantification, and network, proteogenomics, and protein turnover analysis.

# Release notes (Version 2.0.0)
In this version 
1. We assume that the overall amount of proteins in mice is unchanged over time as the mice used in this study are adults. 
2. The Lys used in protein synthesis originate from the food intake, with the rest recycled through protein degradation. 
3. Concentration of soluble Lys is conserved; the rate of free Lys absorbed from food is assumed to be equal to the rate excreted. 
4. Three different settings, which represents simple to comprehensive inputs 

# Software Requirements
The program was written in MATLAB language. It should run on every system with MATLAB R2014 and above.

- MATLAB toolbox needed: 

We recommend using MATLAB R2014 (The MathWorks, Inc., Natick, Massachusetts, United States) or the above version on Linux, Mac or Windows. Our program is mainly using the Global Optimization toolbox other than basic toolboxes.

# Hardware Requirements
The program can be run on either Linux or windows. A memory size of 4GB is optimum, and higher memory may require depending on the size of the data.
The current program has been successfully tested with MATLAB R2019 version on the following system: 16 GB memory 3.3 GHz CPU processors with 6 cores.

# Installation
Installation of the script is not required. Download all the scripts to any working directory (e.g.,/home/usr/JUMPt). IMPORTANT: All the scripts, including associated modules (associated with the program), should be placed in the same folder. 

# Input File Preparation
A testing dataset with 100 proteins (test_data.xlsx) is available along with the scripts for evaluation purposes. Like the testing dataset, the user needs to prepare the input data file with the information below.
1.	pSILAC proteins (mandatory)
2.	pSILAC free (unbolund) Lys (optional)
3.	Free Lys concentration (optional)
4.	Lys concentration in individual proteins (optional)

In the current version, the absolute concentration of each individual protein was calculated based on the commonly used APEX method based on theoretical peptide identification probability and MS/MS spectral counts in a deep brain proteome dataset covering more than 14,000 unique mouse proteins. Each protein-bound Lys concentration was then calculated according to the APEX output, Lys residue number in each protein, and the total measured protein-bound Lys concentration (i.e., 41300 microM). A list of 14000 brain proteins and their concentrations can be found in the file "Brain_proteinConcentrations.xlsx".  

# Run JUMPt program (Demo data set)
JUMPt requires a parameter file (JUMPt.parms). The user needed to specify the following parameters in the 'JUMPt.params' file.
1. JUMPt setting for calculating the half-lives. 
2. Input and output file names (along with the exact path) 
3. Total Absolute Lys concentration from the specific tissue (optional)

Once the parameter file is ready, run "PT_main.m" in MATLAB.

Non-linear fitting of proteins and Lys data using ODE is computationally expensive, especially when the protein data is huge (e.g. > 1000 proteins).  To reduce the computational complexity, we divide the data into sets of 100 proteins. The program performs multistep optimization to find the optimal degradation rates (turnover rates or Half-lives). 
The final output with protein half-lives (in days) and their confidence intervals were saved in the output file mentioned in the params file. The results can be used to understand the turnover rate of the protein of interest.

# Maintainers
To submit bug reports and feature suggestions, please contact:

Surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

# Acknowledgment
We gratefully acknowledge St. Jude Children's Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite.
