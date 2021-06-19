%% program to study the protein turnover kinetics; version 2.0.1 ( last upadted on 05-21-2021)

close all
clearvars 
%% 
%read parameter file
fid     = fopen('JUMPt.params');
param   = textscan(fid,'%s','delimiter','\n');
cellfun(@eval,param{1});
fclose(fid);

params.setting      = setting;
params.input_file   = input_file;
params.bin_size     = bin_size;
params.opti_algo    = optimization_algorithm;
params.purity       = purity_of_SILAC_food;

%Read the input data from input file and distribute the protein different bins 
data = binning(params);
fprintf('\n *******  Completed reding input file; now going to fit protein data and calculate half-lives *******\n\n')

% Calculate the half-lives with different settings
calc_half_lives(data,params); %Fitting/optimizing proteins to calculate the half-lives


fprintf('\n *******  Completed exporting half-lives to out_file ******* \n *******  JUMPt program is complete *******\n\n')
