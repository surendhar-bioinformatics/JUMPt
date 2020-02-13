# Contents of this file
- Introduction
- Release notes
- Software Requirements
- Hardware Requirements
- Installation
- Major steps in the program (Demo data set) 
- Maintainers
- Acknowledgements

# Introduction
JUMPt (JUMP-turnover) software is for determining the protein turnover rates in metabolically labeled animals using mass spectrometry (MS) data. JUMPt uses novel differential equation-based mathematical model to determine the reliable and accurate protein turnover rates. The proposed method determines the half-life of individual proteins by fitting the dynamic data of unlabeled free Lys and protein-bound Lys from individual proteins simultaneously.

JUMPt is part of JUMP Software Suite (shortly JUMP), which is an integrative omics data processing and analysis tool, including protein/peptide database creation, database search, identification filtering, quantification, network analysis, proteogenomics and protein turnover analysis.

# Release notes (Version 1)
In this version 
1. We assume that the overall amount of proteins in mice is unchanged over time as the mice used in this study are adult. 
2. The Lys used in protein synthesis originate from the food intake, with the rest recycled through protein degradation. 
3. oncentration of soluble Lys is conserved, the Lys that is absorbed from the food excreted by degradation pathways. 

# Software Requirements
The program is written in Python3 and MATLAB language. It should run on every system with Python3 with required modules loaded and MATLAB R2017 and above. The minimum required Perl version should be Perl 5.8 or better.

- Python modules are needed:

The following python modules are required. For new user, we encourage to use Anacoda, a free and open-source distribution of the Python.
numpy
scipy
pandas
matplotlib
subprocess 
To install Anaconda distribution, please refer to: https://www.anaconda.com/

- MATLAB toolbox needed: 

We recommend using MATLAB R2014 (The MathWorks, Inc., Natick, Massachusetts, United States) or above version on Linux, Mac or Windows. Our program is mainly using the following toolbox other than basic toolboxes.
Global Optimization toolbox

# Hardware Requirements
The program can be run on either Linux or windows. Memory size of 8GB is an optimum and higher memory may require depending on size of the data.
The current program has been successfully tested on the following system: 16 GB memory 3.3 GHz CPU processors with 6 cores.

# Installation
Installation of the script is not required. Download all the scripts to any working directory (e.g. /home/usr/JUMPt). IMPORTANT: All the scripts including associated modules (associated with the program) should be placed in the same folder. Once the scripts are saved, open and run "main.py" in python console.

# Detailed method of protein turnover calculation:

In order to determine protein degradation rates, we fit the experimental data of all the proteins and Lys using the matrix exponential function. The ODE for Lys (Eq. 11) and individual proteins (Eq. 9) arranged into matrix form  <img src="https://latex.codecogs.com/svg.latex?\Large&space;\frac{d}{dt}\mathbf{\theta}\left(t\right)=\mathbf{G\theta}(t)" title="S1" />
where the matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{G}" /> is unknown with degradation parameters ( <img src="https://render.githubusercontent.com/render/math?math=\gamma_a , \gamma_i \, and \, \gamma_U"> ) to be determined. We rearrange the equation (11) and (9) as follows:


<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{A_L}}{dt}=\gamma_a\left(\theta_{F_L}-\theta_{A_L}\right)\&space;&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\left(\theta_{P_i}-\theta_{A_L}\right)}&plus;\gamma_U\frac{\left[\eta_UP_U\right]}{\left[A\right]}\left(\theta_U-\theta_{A_L}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{A_L}}{dt}=\gamma_a\left(\theta_{F_L}-\theta_{A_L}\right)\&space;&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\left(\theta_{P_i}-\theta_{A_L}\right)}&plus;\gamma_U\frac{\left[\eta_UP_U\right]}{\left[A\right]}\left(\theta_U-\theta_{A_L}\right)"  /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{A_L}}{dt}=\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\theta_{P_i}}-\left(\gamma_a&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}}&plus;\frac{\left[\eta_iP_i\right]}{\left[A\right]}\right)\theta_{A_L}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{A_L}}{dt}=\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\theta_{P_i}}-\left(\gamma_a&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}}&plus;\frac{\left[\eta_iP_i\right]}{\left[A\right]}\right)\theta_{A_L}"  /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{P_i}}{dt}=\gamma_i\left(\theta_{A_L}-\theta_{P_i}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{P_i}}{dt}=\gamma_i\left(\theta_{A_L}-\theta_{P_i}\right)" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{P_i}}{dt}=\gamma_{i\theta_{A_L}}-\gamma_{i\theta_{P_i}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{P_i}}{dt}=\gamma_{i\theta_{A_L}}-\gamma_{i\theta_{P_i}}"  /></a>

We construct our optimization problem with 
- the function <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}:\mathbb{R}\mapsto\mathbb{R}^{n&plus;3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}:\mathbb{R}\mapsto\mathbb{R}^{n&plus;3}" title="\mathbf{\theta}:\mathbb{R}\mapsto\mathbb{R}^{n+3}" /></a> which is built out of the n functions <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{P_{iL}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{P_{iL}}" title="\theta_{P_{iL}}" /></a> along with <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{A_L}&space;and&space;\theta_F" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{A_L}&space;and&space;\theta_F" title="\theta_{A_L} and \theta_F" /></a>.  Without loss of generality, we set <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}_i=&space;\theta_{P_{iL}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}_i=&space;\theta_{P_{iL}}" /></a> for <img src="https://latex.codecogs.com/gif.latex?1\le&space;i\le&space;n" title="1\le i\le n" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}_{n&plus;1}=&space;\theta_{A_L}"  /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}_{n&plus;2}=&space;\theta_F,&space;\mathbf{\theta}_{n&plus;3}=&space;1"  />.  
The n+3th entry is meant to represent the lumped unidentified proteome (U). 

- The vector <img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}\in\mathbb{R}^{n&plus;2}" />, defined similarly to <img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}:&space;\mathbf{\gamma}_i=&space;\gamma_i" /> for <img src="https://latex.codecogs.com/gif.latex?1\le&space;i\le&space;n" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_{n&plus;\mathbf{1}}=\gamma_a" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_{n&plus;\mathbf{2}}=\gamma_U" />. 
 
- The mapping <img src="https://latex.codecogs.com/gif.latex?G:\mathbb{R}^{n&plus;2}\mapsto\&space;\mathbb{R}^{n&plus;3\&space;\times\&space;n&plus;3}" /> which constructs a matrix from a vector <img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}" /> as, 

<a href="https://www.codecogs.com/eqnedit.php?latex={\mathbf{G}\left(\mathbf{\gamma}\right)}_{ij}=\left\{\begin{matrix}&space;{-\gamma}_i&space;&&space;if\&space;i=j\le&space;n\&space;or\&space;i=j=n&plus;3\\&space;-\gamma_i&space;&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}}&plus;\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;if\&space;i=j=n&plus;1\\&space;\gamma_i&space;&&space;if\&space;\left(i\le&space;n\&space;and\&space;j=n&plus;1\right)\&space;or&space;\\&space;&&space;\left(i=n&plus;3\&space;and\&space;j=n&plus;1\right)\\&space;\frac{\left[\eta_i&space;P_i\right]}{\left[A\right]}\gamma_i&space;&&space;if\&space;i=n&plus;1\&space;and\&space;j\le&space;n\\&space;\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;if\&space;i=n&plus;1\&space;and\&space;j=n&plus;3\\&space;\gamma_a&space;&&space;if\&space;i=n&plus;1\&space;and\&space;j=n&plus;2\\&space;0&space;&&space;otherwise&space;\end{matrix}\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\mathbf{G}\left(\mathbf{\gamma}\right)}_{ij}=\left\{\begin{matrix}&space;{-\gamma}_i&space;&&space;if\&space;i=j\le&space;n\&space;or\&space;i=j=n&plus;3\\&space;-\gamma_i&space;&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}}&plus;\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;if\&space;i=j=n&plus;1\\&space;\gamma_i&space;&&space;if\&space;\left(i\le&space;n\&space;and\&space;j=n&plus;1\right)\&space;or&space;\\&space;&&space;\left(i=n&plus;3\&space;and\&space;j=n&plus;1\right)\\&space;\frac{\left[\eta_i&space;P_i\right]}{\left[A\right]}\gamma_i&space;&&space;if\&space;i=n&plus;1\&space;and\&space;j\le&space;n\\&space;\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;if\&space;i=n&plus;1\&space;and\&space;j=n&plus;3\\&space;\gamma_a&space;&&space;if\&space;i=n&plus;1\&space;and\&space;j=n&plus;2\\&space;0&space;&&space;otherwise&space;\end{matrix}\right."  /></a>

When the argument <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}" title="\mathbf{\gamma}" /></a> is unambiguous, we simply write \mathbf{G} instead of \mathbf{G}\left(\mathbf{\gamma}\right).
Now we have \mathbf{G}\ \mathbf{\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}(t), and therefore \mathbf{\theta}\left(t\right)=e^{t\mathbf{G}}\mathbf{\theta}(0), where e^\mathbf{G} is the matrix exponential operator of a matrix \mathbf{G}, and (equation (9) and (11) in the main text) have been aggregated into matrix form.  We formalize our minimization problem as

# Major steps in the program (Demo data set)

A testing dataset (test_data.xlsx) is available along with the scripts for evaluation purposes. Start the "main.py" script in Python console. and follow the instruction. Input your data file (preferably in .xlsx) along with full path. During the execution of program, expected output is saved in "Results" folder which is created in the working directory. Run time for demo on a desktop computer is about 2-10 min depending on the system configuration.

The program performs multistep optimization to find the optimal degradation rates (turnover rates or Half-lives). The final output with protein degradation rates and half-lives were saved in the "Results" folder and temporary output of individual optimization and derivative are also saved in Results folder and will be deleted at the end of program. Three major stages of parameter optimization performed in the program are described below.

- Stage-1: Obtaining local Ansatz for protein degradation rates:

Non-linear fitting of proteins and Lys data using ODE is computationally expensive especially when the protein data is huge (e.g. > 1000). And it is very imprint to provide a best initial guess to global optimization problem to estimate good parameters. In order to get a good initial guess for the protein degradation rates, we fit each protein and Lys using same ODE as we described in the main method. Here we also included rest of the proteome as un-identified protein.
In this step we fit each protein with least residual as possible and calculate the derivates of simulated protein dynamics. The results of initial protein fitting are saved in "Indi_fit_res.xlsx" which is further used in Stage-2 to Globally fit all the proteins and Lys dynamics and obtain degradation rates.

- Stage-2: Global optimization for protein degradation rate using derivatives:

In this stage, we globally fit the derivatives of each individual proteins and Lys data obtained in the previous step. The global fitting of proteins and Lys derivatives will provide us the degradation rates (half-lives) of each protein and free Lys. The results will be saved in "Glob_fit_res.xlsx". One can use these results to understand the turnover rate of protein of interest.

- Stage-3 (Optional): Global optimization of protein and degradation rate:

Global optimization algorithms in MATLAB are efficient to handle large system, but computationally very expensive. If user is interested to further Globally optimize the parameters in MATLAB, they can continue with this stem. In our experimented this step shows a slightly decrease in the over residual error. User can choose whether to continue with this step or not at this stage.


# Maintainers
To submit bug reports and feature suggestions, please contact:

Surendhar Reddy Chepyala (surendharreddy.chepyala@stjude.org) and Junmin Peng (junmin.peng@stjude.org)

# Acknowledgement
We gratefully acknowledge St. Jude Children’s Research Hospital, ALSAC (American Lebanese Syrian Associated Charities) and National Institute of Health for supporting the development of JUMP Software Suite. 
