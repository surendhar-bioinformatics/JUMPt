%% program to study the protein turnover kinetics; version 2.0.1 ( last upadted on 05-21-2021)

close all
clearvars 
%% 
%read parameter file
fid = fopen('JUMPt.params');
param = textscan(fid,'%s','delimiter','\n');
cellfun(@eval,param{1});
fclose(fid);
%% read the data
All_data    = readtable(input_file);
mask        = startsWith(All_data.Properties.VariableNames, 'data');
mask3       = startsWith(All_data.Properties.VariableNames, 'Concentration');
mask2       = startsWith(All_data.Name, 'Lys');
data.t      = table2array(All_data(1,mask));  
data.t_long = linspace(0, 32,321);

if any(mask2)
    data.Lys        = table2array(All_data(mask2,mask));% Lys data
    data.protInfo   = (All_data(3:end,1:find(mask,1)-1)); % protein IDs
    data.ProtConc   = table2array(All_data(3:end,mask3))';% Protein concentration
    SILAC_data      = table2cell(All_data(3:end,mask));% protein data
    SILAC_data_NaN  = cellfun(@(x) ~isa(x,'double'),SILAC_data); % check for possible 0x0 char reads
    data.LysConc  = table2array(All_data(2,mask3));
    SILAC_data(SILAC_data_NaN)   = {NaN};% set those to NaN
    data.SILAC_data              = cell2mat(SILAC_data)'; % heavy/light ratios,
else
    data.protInfo  = (All_data(2:end,1:find(mask,1)-1));
    data.ProtConc  = table2array(All_data(2:end,mask3))'; 
    SILAC_data     = table2cell(All_data(2:end,mask));
    SILAC_data_NaN = cellfun(@(x) ~isa(x,'double'),SILAC_data); 
    SILAC_data(SILAC_data_NaN)   = {NaN};
    data.SILAC_data              = cell2mat(SILAC_data)'; 
end
parameters.setting1 = Setting_1;
parameters.setting2 = Setting_2;
parameters.setting3 = Setting_3;
parameters.tot_Lys  = Total_Lys_concentration_in_tissue;
parameters.out_file = out_file;
parameters.prot_per_fit     = Num_prot_in_each_fitting;
parameters.Num_prot_to_fit  = size(data.SILAC_data,2);

%calling The  main function
Prot_turnover_v1(data,parameters)

fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')


