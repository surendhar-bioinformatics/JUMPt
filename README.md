# Contents of this file
- Introduction
- Release notes
- Software Requirements
- Hardware Requirements
- Installation
- Major steps in the program  
- Maintainers
- Acknoledgements

# Introduction
JUMPt (JUMP-turnover) software is for detrmining the protein turnover rates in metabolically labeled animals using mass spectrometry(MS) data. JUMPt uses novel differential equation-based mathematical model to determine the reliable and accurate protein turnover rates. The proposed method determines the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously. 

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

The program perform multistep optimization to find the optimal degredation rates (turnover rates or Half-lives). The final output with protein degredation rates and halflives were saved in the "Results" folder and temperory output of individual optimization and derivatie are also saved in Results folder and will be deleted at the end of program. Three major stages of parameter optimization performed in the program are described below.

- Stage-1: Obtaining local Ansatz for protein degrdation rates:

Non-linear fitting of proteins and Lys data using ODE is computationally expensive especially when the protein data is huge ( eg. > 1000). And it is very impront to provide a best inial guess to global opimization problem  to estimate a good paramters. In order to get a good inial guess for the protein degredation rates, we fit the each protein and Lys using same ODE as we described in the main method. Here we also included rest of the proteome as un-identified protein. 

In this step we fit each protein with least residual as possible and calculate the derivates of simulated protein dynamics. The results of inditial protein fitting are saved in "Indi_fit_res.xlsx" which is further used in Stage-2 to Globally fit all the proteins and Lys dynamics and obtain degredation rates. 

- Stage-2: Global optimization for protein degredation rate using derivatives:

In this stage, we globally fit the derivatives of each individual proteins and Lys data obtained in the previous step. The global fitting of proteins and Lys derivatives will provide us the degredation rates (half-lifes) of each protein and free Lys. The results will be saved in "Glob_fit_res.xlsx". One can use these results to understnad the turnover rate of protein of intrest.

- Stage-3 (Optional): Global optimization of protein and degredation rate:
Global optimization algorithms in MATLAB are efficiant to handle large system, but computationally very expensive. If user is intrested to further Globally optimize the parameters in MATLAB, they can ccontinue with this stem. In our experimance this step shows a slghly decrese in the over residual error. User can choose whether to continue with this step or not at this stage.


# Maintainers
To submit bug reports and feature suggestions, please contact:

surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

# Acknowledgement
We gratefully acknowledge St. Jude Children’s Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite. 
