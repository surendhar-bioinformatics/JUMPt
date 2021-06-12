%% program to study the protein turnover kinetics; version 2.0.1 ( last upadted on 05-21-2021)

close all
clearvars 
%% 
%read parameter file
fid = fopen('JUMPt.params');
param = textscan(fid,'%s','delimiter','\n');
cellfun(@eval,param{1});
fclose(fid);

parameters.setting1 = Setting_1;
parameters.setting2 = Setting_2;
parameters.setting3 = Setting_3;
parameters.tot_Lys  = Total_Lys_concentration_in_tissue;
parameters.input_file = input_file;
parameters.out_file   = out_file;
parameters.prot_per_fit     = Num_prot_in_each_fitting;

%Split the proteins into bin
data = Binning(parameters);

%% Calculate the half-lives with different settings

Glob_fit_Prot = [];
if parameters.setting3 == 1
	Glob_fit_Prot_setting3 = gama_Prot(data,parameters,[], 3); %Fitting/optimizing proteins to calculate the half-lives
end

if parameters.setting2 == 1
    OptiResStep1 = gama_Lys(data,parameters,2);% Fitting/optimizing proteins with all time points to get the Lys degredation rate 
    if (data.se2(end) ~= data.se(end)) || (length(data.se) ~= 2)
        Glob_fit_Prot_setting2 = gama_Prot(data,parameters,OptiResStep1, 2); % Fitting/optimizing proteins to calculate the half-lives
    end
end

if parameters.setting1 == 1
    OptiResStep1 = gama_Lys(data,parameters,1); 
    if (data.se2(end) ~= data.se(end)) || (length(data.se) ~= 2)
        Glob_fit_Prot_setting1 = gama_Prot(data,parameters,OptiResStep1,1);
    end
end

fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')


