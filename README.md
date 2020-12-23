# Contents of this file
- Introduction
- Release notes
- Software Requirements
- Hardware Requirements
- Installation
- Input File Preparation
- Run JUMPt program (test data provided) 
- Maintainers
- Acknowledgements

# Introduction
JUMPt (JUMP-turnover) software is for determining the protein turnover rates in metabolically labeled animals using mass spectrometry (MS) data. JUMPt uses novel differential equation-based mathematical model to determine the reliable and accurate protein turnover rates. The proposed method determines the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously.

JUMPt is part of JUMP Software Suite (shortly JUMP), which is an integrative omics data processing and analysis tool, including protein/peptide database creation, database search, identification filtering, quantification, network analysis, proteogenomics and protein turnover analysis.

# Release notes (Version 1.0.0)
In this version 
1. We assume that the overall amount of proteins in mice is unchanged over time as the mice used in this study are adult. 
2. The Lys used in protein synthesis originate from the food intake, with the rest recycled through protein degradation. 
3. Concentration of soluble Lys is conserved, the rate of free Lys absorbed from food is assumed to be equal to the rate excreted. 

# Software Requirements
The program is written in MATLAB language. It should run on every system with MATLAB R2014 and above.

- MATLAB toolbox needed: 

We recommend using MATLAB R2014 (The MathWorks, Inc., Natick, Massachusetts, United States) or above version on Linux, Mac or Windows. Our program is mainly using the Global Optimization toolbox other than basic toolboxes.

# Hardware Requirements
The program can be run on either Linux or windows. Memory size of 4GB is an optimum and higher memory may require depending on size of the data.
The current program has been successfully tested on the following system: 16 GB memory 3.3 GHz CPU processors with 6 cores.

# Input File Preparation
Installation of the script is not required. Download all the scripts to any working directory (e.g. /home/usr/JUMPt). IMPORTANT: All the scripts including associated modules (associated with the program) should be placed in the same folder. Once the scripts are saved, open and run "PT_main.m" in MATLAB.
JUMPt program requires input data file (preferably in .xlsx) with following information.
1.	SILAC ratio for proteins and free Lys
2.	Free Lys and total Lys concentrations
3.	Lys concentration in in individual proteins which is calculated based on the proteins concentration multiplied by number of Lys molecule in it.

# Run JUMPt program (Demo data set)
A testing dataset (test_data.xlsx) is available along with the scripts for evaluation purposes. Start the " PT_main.m " script in MATLAB. and follow the instruction. The work flow of the JUMPt software is describe din Figure 1. 

Non-linear fitting of proteins and Lys data using ODE is computationally expensive especially when the protein data is huge (e.g. > 1000).  To reduce the complexity of we do devide the data into sets of 100 proteins. For each set the program performs multistep optimization to find the optimal degradation rates (turnover rates or Half-lives). 
The final output with protein degradation rates and half-lives were saved in the same working directory. One can use these results to understand the turnover rate of protein of interest.

Global optimization of protein and degradation rate:


# Maintainers
To submit bug reports and feature suggestions, please contact:

Surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

# Acknowledgement
We gratefully acknowledge St. Jude Children’s Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite. 
