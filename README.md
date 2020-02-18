# Contents of this file
- Introduction
- Release notes
- Software Requirements
- Hardware Requirements
- Installation
- Details of protein turnover calculation method
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

# Details of protein turnover calculation method

In order to determine protein degradation rates, we fit the experimental data of all the proteins and Lys using the matrix exponential function. The ODE for Lys and individual proteins (equation(11) and (9)) can be rewritten in matrix form  <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\frac{d}{dt}\mathbf{\theta}(t)&space;=&space;\mathbf{G}\&space;\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\frac{d}{dt}\mathbf{\theta}(t)&space;=&space;\mathbf{G}\&space;\mathbf{\theta}\left(t\right)" /></a>
where the matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{G}" /> has degradation parameters ( <img src="https://render.githubusercontent.com/render/math?math=\gamma_a , \gamma_i \, and \, \gamma_U"> ) to be determined. We rearrange the equation (9) and (11) in two different forms as follows:


<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{A_L}}{dt}=\gamma_a\left(\theta_{F_L}-\theta_{A_L}\right)\&space;&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\left(\theta_{P_i_L}-\theta_{A_L}\right)}&plus;\gamma_U\frac{\left[\eta_UP_U\right]}{\left[A\right]}\left(\theta_{U_L}-\theta_{A_L}\right)," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{A_L}}{dt}=\gamma_a\left(\theta_{F_L}-\theta_{A_L}\right)\&space;&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\left(\theta_{P_i_L}-\theta_{A_L}\right)}&plus;\gamma_U\frac{\left[\eta_UP_U\right]}{\left[A\right]}\left(\theta_{U_L}-\theta_{A_L}\right)," /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\small&space;\frac{{d\theta}_{A_L}}{dt}=\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\theta_{P_i_L}}-\left(\gamma_a&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}}&plus;\frac{\left[\eta_iP_i\right]}{\left[A\right]}\right)\theta_{A_L}&plus;\gamma_a\theta_{F_L}&plus;\gamma _U\frac{\left[\eta_UP_U\right]}{\left[A\right]}\theta_{U_L}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\small&space;\frac{{d\theta}_{A_L}}{dt}=\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}\theta_{P_i_L}}-\left(\gamma_a&plus;\sum_{i=1}^{n}{\gamma_i\frac{\left[\eta_iP_i\right]}{\left[A\right]}}&plus;\frac{\left[\eta_iP_i\right]}{\left[A\right]}\right)\theta_{A_L}&plus;\gamma_a\theta_{F_L}&plus;\gamma_U\frac{\left[\eta_UP_U\right]}{\left[A\right]}\theta_{U_L}," /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{P_i_L}}{dt}=\gamma_{i\theta_{A_L}}-\gamma_{i\theta_{P_i_L}}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{P_i_L}}{dt}=\gamma_{i\theta_{A_L}}-\gamma_{i\theta_{P_i_L}},"  /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{{d\theta}_{P_i_L}}{dt}=\gamma_{i}\theta_{A_L}-\gamma_{i}\theta_{P_{iL}}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{{d\theta}_{P_i_L}}{dt}=\gamma_{i}\theta_{A_L}-\gamma_{i}\theta_{P_{iL}}." /></a>

We construct our optimization problem with 
- the function <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{\theta}:\mathbb{R_&plus;}\mapsto\mathbb{R}_&plus;^{n&plus;3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbf{\theta}:\mathbb{R_&plus;}\mapsto\mathbb{R}_&plus;^{n&plus;3}" /></a> 
which is built out of the n functions <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{P_{iL}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{P_{iL}}" title="\theta_{P_{iL}}" /></a> along with <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{A_L},&space;\&space;\theta_F_L&space;\&space;and&space;\&space;\theta_{U_L}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{A_L},&space;\&space;\theta_F_L&space;\&space;and&space;\&space;\theta_{U_L}" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbb{R_&plus;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbb{R_&plus;}" /></a>  is the set of nonnegative real numbers.  Without loss of generality, we set <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}_i=&space;\theta_{P_{iL}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}_i=&space;\theta_{P_{iL}}" /></a> for <img src="https://latex.codecogs.com/gif.latex?1\le&space;i\le&space;n" title="1\le i\le n" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}_{n&plus;1}=&space;\theta_{A_L}"  /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}_{n&plus;2}=&space;\theta_F_L,&space;\mathbf{\theta}_{n&plus;3}=&space;\theta_U_L"  />.   The n+3th entry is meant to represent the lumped unidentified proteome (U). 

- The vector <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{g}\in\mathbb{R}_&plus;^{n&plus;2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbf{g}\in\mathbb{R}_&plus;^{n&plus;2}" /></a> defined similarly to <img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}:&space;\mathbf{g}_i=&space;\gamma_i" /> for <img src="https://latex.codecogs.com/gif.latex?1\le&space;i\le&space;n" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{g}_{n&plus;\mathbf{1}}=\gamma_a" /> and <img src="https://latex.codecogs.com/gif.latex?\mathbf{g}_{n&plus;\mathbf{2}}=\gamma_U" />. 
 
- The mapping <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\mathbf{G}:\mathbb{R}_&plus;^{n&plus;2}\mapsto&space;\mathbb{R}_&plus;^{(n&plus;2)\times(n&plus;3)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\mathbf{G}:\mathbb{R}_&plus;^{n&plus;2}\mapsto&space;\mathbb{R}_&plus;^{(n&plus;2)\times(n&plus;3)}" /></a> which constructs a matrix from a vector <img src="https://latex.codecogs.com/gif.latex?\mathbf{g}" /> as, 

<a href="https://www.codecogs.com/eqnedit.php?latex={\mathbf{G}\left(\mathbf{g}\right)}_{ij}=\left\{\begin{matrix}&space;\mathbf{-g}_i&space;&&space;\textup{if}\&space;i=j\le&space;n\&space;\textup{or}&space;\&space;i=j=n&plus;3&space;\\&space;\mathbf{-g_i}&space;-\sum_{i=1}^{n}{\frac{\left[\eta_iP_i\right]}{\left[A\right]}}-\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;\textup{if}\&space;i=j=n&plus;1\\&space;\mathbf{g}_i&space;&&space;\textup{if}\&space;(i\le&space;n\&space;and\&space;j=n&plus;1)&space;\&space;\textup{or}&space;\\&space;\&space;&&space;(i=n&plus;3\&space;\textup{and}\&space;j=n&plus;1)\\&space;\frac{\left[\eta_j&space;P_j\right]}{\left[A\right]}\mathbf{g}_i&space;&&space;\textup{if}\&space;i=n&plus;1\&space;\textup{and}\&space;j\le&space;n\\&space;\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;\textup{if}\&space;i=n&plus;1\&space;\textup{and}\&space;j=n&plus;3\\&space;\mathbf{g}_i&space;&&space;\textup{if}\&space;i=n&plus;1\&space;\textup{and}\&space;j=n&plus;2\\&space;0&space;&&space;\textup{otherwise}&space;\end{matrix}\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\mathbf{G}\left(\mathbf{g}\right)}_{ij}=\left\{\begin{matrix}&space;\mathbf{-g}_i&space;&&space;\textup{if}\&space;i=j\le&space;n\&space;\textup{or}&space;\&space;i=j=n&plus;3&space;\\&space;\mathbf{-g_i}&space;-\sum_{i=1}^{n}{\frac{\left[\eta_iP_i\right]}{\left[A\right]}}-\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;\textup{if}\&space;i=j=n&plus;1\\&space;\mathbf{g}_i&space;&&space;\textup{if}\&space;(i\le&space;n\&space;and\&space;j=n&plus;1)&space;\&space;\textup{or}&space;\\&space;\&space;&&space;(i=n&plus;3\&space;\textup{and}\&space;j=n&plus;1)\\&space;\frac{\left[\eta_j&space;P_j\right]}{\left[A\right]}\mathbf{g}_i&space;&&space;\textup{if}\&space;i=n&plus;1\&space;\textup{and}\&space;j\le&space;n\\&space;\frac{\left[\eta_U&space;P_U\right]}{\left[A\right]}&space;&&space;\textup{if}\&space;i=n&plus;1\&space;\textup{and}\&space;j=n&plus;3\\&space;\mathbf{g}_i&space;&&space;\textup{if}\&space;i=n&plus;1\&space;\textup{and}\&space;j=n&plus;2\\&space;0&space;&&space;\textup{otherwise}&space;\end{matrix}\right." /></a>

When the argument <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{g}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{g}" /></a> is unambiguous, we simply write <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G}" /></a> instead of <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G}\left(\mathbf{g}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G}\left(\mathbf{g}\right)" /></a>.

Now we have <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G}\&space;\mathbf{\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}(t)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G}\&space;\mathbf{\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}(t)" /></a>, and therefore <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}\left(t\right)=e^{t\mathbf{G}}\mathbf{\theta}(0)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}\left(t\right)=e^{t\mathbf{G}}\mathbf{\theta}(0)" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=e^\mathbf{G}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?e^\mathbf{G}" /></a> is the matrix exponential operator of a matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G}" /></a>, and (equation (9) and (11) in the main text) have been aggregated into matrix form.  We formalize our minimization problem as

<<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{g}_\ast={\rm&space;argmin,\&space;}\below{\mathbf{g}\in\mathbf{\mathbb{R}}^{n&plus;2}}&space;\sum_{i=1}^{m}\left&space;\|&space;e^{t_i\mathbf{G}}\mathbf{\theta}\left(0\right)-\mathbf{\theta}\left(t_i\right)\right&space;\|_{W}^{2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{g}_\ast={\rm&space;argmin,\&space;}\below{\mathbf{g}\in\mathbf{\mathbb{R}}^{n&plus;2}}&space;\sum_{i=1}^{m}\left&space;\|&space;e^{t_i\mathbf{G}}\mathbf{\theta}\left(0\right)-\mathbf{\theta}\left(t_i\right)\right&space;\|_{W}^{2}" /></a>

where the <a href="https://www.codecogs.com/eqnedit.php?latex=t_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t_i" /></a> are the time points for which we have protein data and <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}\left(0\right)=1," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}\left(0\right)=1" /></a> is our initial condition.  The matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{W}\in\mathbb{R}^{(n&plus;3),\times&space;(n&plus;3)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{W}\in\mathbb{R}^{(n&plus;3),\times&space;(n&plus;3)}," /></a>, defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{W}_{ij}=\left\{\begin{matrix}1&\textup{if}\&space;i=j\le&space;n&plus;2\\0&\textup{otherwise}\\\end{matrix}\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{W}_{ij}=\left\{\begin{matrix}1&\textup{if}\&space;i=j\le&space;n&plus;2\\0&\textup{otherwise}\\\end{matrix}\right." /></a>

which exists solely to exclude the dimensions of <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}\left(t\right)" /></a> that are constructed from <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{A_F}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{A_F}" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_U_L" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_U_L"  /></a> from calculation of error.

Though the optimal degradation parameters can be found with a gradient-based search of the n+2 dimensional positives <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{R}^{n&plus;2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{R}^{n&plus;2}"  /></a>.  Due to the mapping <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbf{G}" /> and the matrix exponential, a closed-form expression for the gradient <a href="https://www.codecogs.com/eqnedit.php?latex=\nabla\mathbf{\gamma}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nabla\mathbf{\gamma}" /></a> is not available, and thus we avoid gradient-based search to find the optimal <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_\ast" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_\ast" /></a>.  Instead, we propose derivative fitting approach and the details of the algorithm are given below:


1. Fit <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}\left(t\right)"  /></a> with a surrogate model <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vartheta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vartheta}\left(t\right)" /></a> so that <a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;\|&space;\mathbf{\theta}\left(t_i\right)-\mathbf{\vartheta}\left(t_i\right)&space;\right&space;\|_{W}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left&space;\|&space;\mathbf{\theta}\left(t_i\right)-\mathbf{\vartheta}\left(t_i\right)&space;\right&space;\|_{W}" /></a> is minimized over all <a href="https://www.codecogs.com/eqnedit.php?latex=t_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t_i" /></a>  and both  <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vartheta}\left(t\right)&space;\&space;and&space;\&space;\frac{d}{dt}\mathbf{\vartheta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vartheta}\left(t\right)&space;\&space;and&space;\&space;\frac{d}{dt}\mathbf{\vartheta}\left(t\right)" /></a> are cheaply available for any  <a href="https://www.codecogs.com/eqnedit.php?latex=t\in\mathbf{P}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?t\in\mathbf{P}." /></a>
 We fit each protein <a href="https://www.codecogs.com/eqnedit.php?latex=(\mathbf{\theta}_{\mathbit{P}_{\mathbit{iL}}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(\mathbf{\theta}_{\mathbit{P}_{\mathbit{iL}}})"  /></a> individually together with Lys <a href="https://www.codecogs.com/eqnedit.php?latex=(\mathbf{\theta}_{\mathbit{AL}})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(\mathbf{\theta}_{\mathbit{AL}})" /></a>  data to get both  <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vartheta}\left(t\right)&space;and&space;\frac{d}{dt}\mathbf{\vartheta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vartheta}\left(t\right)&space;and&space;\frac{d}{dt}\mathbf{\vartheta}\left(t\right)"  /></a>.

2. Pick a uniform grid of  k points <a href="https://www.codecogs.com/eqnedit.php?latex=\widehat{t_i}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widehat{t_i}"  /></a> over <a href="https://www.codecogs.com/eqnedit.php?latex=\left[t_1,t_k\right]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left[t_1,t_k\right]" /></a>
	
3. Set j≔1

4. Set the (n+3) x (n+2) matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{M}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{M}"  /></a> as 

<a href="https://www.codecogs.com/eqnedit.php?latex={\mathbf{M}\left(t\right)}_{ij}\left\{\begin{matrix}&space;\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{1}}\left(t\right)-\mathbf{\vartheta}_\mathbit{i}\left(t\right)&space;&&space;if\&space;i=j\le\&space;n\\&space;\left(\mathbf{\vartheta}_\mathbit{i}\left(t\right)-\&space;\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{1}}\left(t\right)\right)&&space;if\&space;i=n&plus;1\&space;and\&space;j\le&space;n\\&space;\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{2}}\left(t\right)-\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{1}}\left(t\right)&&space;if\&space;i=n&plus;1\&space;and\&space;j=n&plus;1&space;\\&space;0&space;&&space;otherwise&space;\\&space;\end{matrix}\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\mathbf{M}\left(t\right)}_{ij}\left\{\begin{matrix}&space;\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{1}}\left(t\right)-\mathbf{\vartheta}_\mathbit{i}\left(t\right)&space;&&space;if\&space;i=j\le\&space;n\\&space;\left(\mathbf{\vartheta}_\mathbit{i}\left(t\right)-\&space;\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{1}}\left(t\right)\right)&&space;if\&space;i=n&plus;1\&space;and\&space;j\le&space;n\\&space;\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{2}}\left(t\right)-\mathbf{\vartheta}_{\mathbit{n}&plus;\mathbf{1}}\left(t\right)&&space;if\&space;i=n&plus;1\&space;and\&space;j=n&plus;1&space;\\&space;0&space;&&space;otherwise&space;\\&space;\end{matrix}\right."  /></a>


5. Let

<a href="https://www.codecogs.com/eqnedit.php?latex=h\left(x\right)=\left\{\begin{matrix}x^2&if\&space;x<0\\0&otherwise\\\end{matrix}\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?h\left(x\right)=\left\{\begin{matrix}x^2&if\&space;x<0\\0&otherwise\\\end{matrix}\right." /></a>

6. set

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_0,u&space;:=&space;{\rm&space;argmin,\&space;}\below{\mathbf{\gamma}\in\mathbb{R}^{n&plus;2};u\in\mathbb{R}^k}=" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_0,u&space;:=&space;{\rm&space;argmin,\&space;}\below{\mathbf{\gamma}\in\mathbb{R}^{n&plus;2};u\in\mathbb{R}^k}=" /></a> <img src="https://render.githubusercontent.com/render/math?math=\sum_{i=1}^{m}\left \| \mathbf{M}\left(\widehat{t_i}\right)\mathbf{\gamma} - \frac{d}{dt}\mathbf{\vartheta}\left(\widehat{t_i}\right) {&plus;} u_i \frac{\left[A\right]}{\left[\eta_U\right]} \right \|_{W}^{2}">

where <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_0,u" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_0,u"  /></a> are solved approximately with an iterative method.



7. Let <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\delta}_{n&plus;3\&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\delta}_{n&plus;3\&space;}" /></a> be the Kronecker delta vector which is zero everywhere except at entry n+3

8. Set the (n+3)rd entry of <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_0" /></a> to zero.  Let <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbit{S}:=&space;\left&space;\{&space;\psi\left(t\right)=e^{\widehat{t_i}\mathbf{G}\left(\mathbf{\gamma}_0&plus;\mathbf{\delta}_{n&plus;3\&space;}g\right)}&space;\mathbf{\theta}\left(0\right):g\geq0&space;\right&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbit{S}:=&space;\left&space;\{&space;\psi\left(t\right)=e^{\widehat{t_i}\mathbf{G}\left(\mathbf{\gamma}_0&plus;\mathbf{\delta}_{n&plus;3\&space;}g\right)}&space;\mathbf{\theta}\left(0\right):g\geq0&space;\right&space;\}" /></a>
be the set of protein dynamics obtained by perturbing <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma_{U}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma_{U}" /></a> and find the <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{U_L}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{U_L}" /></a> that maximizes correlation with <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{u}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{u}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{U_L}:=&space;{\rm&space;\displaystyle&space;argmax,&space;\&space;}\below{\cos{\angle(}\mathbf{u},\&space;\psi)}=&space;{\rm&space;\displaystyle&space;argmax,&space;\&space;}\below{\psi\in\mathbit{S}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{U_L}:=&space;{\rm&space;\displaystyle&space;argmax,&space;\&space;}\below{\cos{\angle(}\mathbf{u},\&space;\psi)}=&space;{\rm&space;\displaystyle&space;argmax,&space;\&space;}\below{\psi\in\mathbit{S}}" /></a> <img src="https://render.githubusercontent.com/render/math?math=\frac{\sum_{i=1}^{k}{u_i\psi\left(\widehat{t_i}\right)}}{(\sum_{i=1}^{k}{u_i^2)(\sum_{i=1}^{k}{\psi\left(\widehat{t_i}\right)^2)}}}" title="\frac{\sum_{i=1}^{k}{u_i\psi\left(\widehat{t_i}\right)}}{(\sum_{i=1}^{k}{u_i^2)(\sum_{i=1}^{k}{\psi\left(\widehat{t_i}\right)^2)}}}">




9. set<a href="https://www.codecogs.com/eqnedit.php?latex={\&space;\mathbf{\gamma}}_0&space;:=&space;\mathbf{\gamma}_0&plus;&space;\mathbf{\delta}_{n&plus;3\&space;}\gamma_{U}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\&space;\mathbf{\gamma}}_0&space;:=&space;\mathbf{\gamma}_0&plus;&space;\mathbf{\delta}_{n&plus;3\&space;}\gamma_{U}" /></a>  where <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma_{U}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma_{U}" /></a> is the perturbation that produced <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_{U_L}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_{U_L}" /></a>  and set 
<a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon_j&space;:=&space;\left&space;\|&space;\mathbf{\theta}\left(t_i\right)-e^{t_i\mathbf{G}\left(\mathbf{\gamma}_0\right)}\mathbf{\theta}\left(0\right)&space;\right&space;\|_{\mathbf{W}}^{2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon_j&space;:=&space;\left&space;\|&space;\mathbf{\theta}\left(t_i\right)-e^{t_i\mathbf{G}\left(\mathbf{\gamma}_0\right)}\mathbf{\theta}\left(0\right)&space;\right&space;\|_{\mathbf{W}}^{2}"  /></a>


10. if j=1 or <a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon_j<\epsilon_{j-1\&space;}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon_j<\epsilon_{j-1\&space;}," /></a> set <a href="https://www.codecogs.com/eqnedit.php?latex=j\&space;:=j&plus;1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?j\&space;:=j&plus;1" /></a>
and goto 6, using <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_0" /></a> as an initial guess for the iterative solver used on line 6.

11. return <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_0" /></a>

**There are some helpful clarifying remarks to be made about Algorithm  we described above**

- We use a grid spacing of 10<sup>-1</sup> for the <a href="https://www.codecogs.com/eqnedit.php?latex=\widehat{t_i}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widehat{t_i}" /></a>  in line 2.
- The matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{M}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{M}\left(t\right)" /></a> on line 4 is a rearrangement of <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}\left(t\right)" /></a>, so it satisfies
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{M}\left(t\right)\mathbf{\gamma}=\frac{d}{dt}\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{M}\left(t\right)\mathbf{\gamma}=\frac{d}{dt}\mathbf{\theta}\left(t\right)" /></a>

- The vector <a href="https://www.codecogs.com/eqnedit.php?latex=u_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_i" /></a> on line 6 is meant to represent a combination of noise and the derivative of the unidentified protein for which no dynamics history is available except for at <a href="https://www.codecogs.com/eqnedit.php?latex=t=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t=0" /></a>.

- Solving the minimization problem on line 6 is done with constrained iterative optimization to enforce <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}\in\mathbb{R}^{n&plus;2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}\in\mathbb{R}^{n&plus;2}" /></a>.  The trust-region method we used occasionally encountered difficulty satisfying all constraints without taking many hundreds of iterations.

- The penalty term <a href="https://www.codecogs.com/eqnedit.php?latex=\tau\left(h\left(u_i\right)\right)^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau\left(h\left(u_i\right)\right)^2" /></a> in the minimization problem on line 6 is intended to penalize <a href="https://www.codecogs.com/eqnedit.php?latex=u\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u\left(t\right)"/></a> that have infeasible dynamics.  
The parameter <a href="https://www.codecogs.com/eqnedit.php?latex=\tau" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau" /></a> is chosen to be the largest possible <a href="https://www.codecogs.com/eqnedit.php?latex=\tau\in\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\tau\in\mathbb{R}" /></a> such that the constraints in the previous remark are not violated. 

- Steps 7--9 find the degredation rates that induces protein dynamics that most closely fit <a href="https://www.codecogs.com/eqnedit.php?latex=u_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_i" /></a>.

The parameters exploits the derivative relationship <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}\left(t\right)" /></a>, where the observed protein dynamics <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\theta}\left(t\right)" /></a> are replaced with a surrogate model <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\vartheta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\vartheta}\left(t\right)" /></a> that is fit to the protein dynamics.  
The matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{M}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{M}\left(t\right)" /></a> simply rearranges <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{G\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}\left(t\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{G\theta}\left(t\right)=\frac{d}{dt}\mathbf{\theta}\left(t\right)"  /></a> so that the vector <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}" /></a>is on the right-hand side of a matrix with known values.  We find the <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{\gamma}_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{\gamma}_0" /></a> that minimizes the L2-norm of the residual over all synthetic points <a href="https://www.codecogs.com/eqnedit.php?latex=\widehat{t_i}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widehat{t_i}" /></a>.

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
