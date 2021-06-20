# Contents
- Introduction
- Release notes
- Software and Hardware Requirements
- Installation
- Input File Preparation
- Run JUMPt program (test data provided) 
- Maintainers
- Acknowledgments
- References

# Introduction
JUMPt (JUMP-turnover) software determines the protein turnover rates in metabolically labeled animals using mass spectrometry (MS) data. JUMPt uses a novel differential equation-based mathematical model to calculate the reliable and accurate protein turnover rates. The proposed method calculates the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously.

JUMPt is part of JUMP Software Suite (shortly JUMP), is an ongoing large software suite developed for the need of mass spectrometry (MS)- based proteomics, metabolomics, and the integration with genomics for network analysis at the level of systems biology. Currently, JUMP can handle protein/peptide database creation, database search, identification filtering, quantification, and network, proteogenomics, and protein turnover analysis.

# Release notes (Version 2.0.0)
In the current version 
1. We assume that the overall amount of proteins in mice is unchanged over time as the mice used in this study are adults. 
2. The Lys used in protein synthesis originate from the food intake, with the rest recycled through protein degradation. 
3. Concentration of soluble Lys is conserved; the rate of free Lys absorbed from food is assumed to be equal to the rate excreted. 
4. Three different settings, which represents simple to comprehensive inputs 

# Software  and Hardware Requirements
The program was written in MATLAB language. The program runs on any Linux, Mac, or Windows computer with MATLAB R2014 (The MathWorks, Inc., Natick, Massachusetts, United States) or the above version. The current JUMPt program has been successfully tested with MATLAB R2019 version on the following system: 16 GB memory 3.3 GHz CPU processors with 6 cores. The program needs more time to complete on the system with fewer core processors in the CPU.
MATLAB toolbox needed: 
- Global Optimization toolbox along with other basic toolboxes

# Installation
Installation of the script is not required. Download all the scripts to any working directory (e.g.,/home/usr/JUMPt). IMPORTANT: All the scripts, including associated modules (associated with the program), should be placed in the same folder. 

# Input File Preparation
A testing dataset with 100 proteins is available for each setting, along with the scripts for evaluation purposes. Similar to the testing dataset, the user needs to prepare the input data file with the information below.

1.	pSILAC proteins (mandatory)
2.	pSILAC free (unbolund) Lys (optional)
3.	Free Lys concentration (optional)
4.	Lys concentration in individual proteins (optional)

In the current version, the absolute concentration of each individual protein was calculated based on the commonly used APEX method based on theoretical peptide identification probability and MS/MS spectral counts in a deep brain proteome dataset covering more than 14,000 unique mouse proteins. Each protein-bound Lys concentration was then calculated according to the APEX output, Lys residue number in each protein, and the total measured protein-bound Lys concentration (i.e., 41300 microM). A list of 14000 brain proteins and their concentrations can be found in the file "Brain_proteinConcentrations.xlsx".  


# Run JUMPt program (Demo data set)

Step-1: Modify the parameter file.
JUMPt requires a parameter file (JUMPt.parms). The user needed to specify the following parameters in the 'JUMPt.params' file.

1.	JUMPt setting 
2.	Input file name (along with the exact path)
3.	Bin size('bin_size') to fit proteins each time 
4.	MATLAB optimization algorithm
5.	Purity of SILAC food 

Step-2: Once the parameter file is ready, run "PT_main.m" in MATLAB.
Non-linear fitting of proteins and Lys data using ODE is computationally expensive, especially when the protein data is huge (e.g.,> 1000 proteins).  To reduce the computational complexity, we divide the proteins into sets with bin sizes between 100 -10. The program finds the optimal degradation rates (turnover rates or Half-lives) by fitting protein data (in setting-1) and free-Lys data (in setting-2 and setting-3).  

An output excel file is generated with a prefix 'results_' to the input file name. The final results with protein half-lives (in days) and their confidence intervals were saved to the output file. In addition, parameters used to calculate the half-lives were saved in the output file. 

# Maintainers

To submit bug reports and feature suggestions, please contact:

Surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

# Acknowledgment

We acknowledge St. Jude Children's Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite.

# References

1.	Chepyala et al.,  JUMPt: Comprehensive protein turnover modeling of in vivo pulse SILAC data by ordinary differential equations. Analytical chemistry (under review)
2.	Wang, X., et al., JUMP: a tag-based database search tool for peptide identification with high sensitivity and accuracy. Molecular & Cellular Proteomics, 2014. 13(12): p. 3663-3673.
3.	Wang, X., et al., JUMPm: A Tool for Large-Scale Identification of Metabolites in Untargeted Metabolomics. Metabolites, 2020. 10(5): p. 190.
4.	Li, Y., et al., JUMPg: an integrative proteogenomics pipeline identifying unannotated proteins in human brain and cancer cells. Journal of proteome research, 2016. 15(7): p. 2309-2320.
5.	Tan, H., et al., Integrative proteomics and phosphoproteomics profiling reveals dynamic signaling networks and bioenergetics pathways underlying T cell activation. Immunity, 2017. 46(3): p. 488-503.
6.	Peng, J., et al., Evaluation of multidimensional chromatography coupled with tandem mass spectrometry (LC/LC− MS/MS) for large-scale protein analysis: the yeast proteome. Journal of proteome research, 2003. 2(1): p. 43-50.
7.	Niu, M., et al., Extensive peptide fractionation and y 1 ion-based interference detection method for enabling accurate quantification by isobaric labeling and mass spectrometry. Analytical chemistry, 2017. 89(5): p. 2956-2963.

