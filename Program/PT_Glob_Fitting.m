%% program to study the protein turnover kinetics. 5/20/2019

close all
clearvars 

%% read the original data
input_file  = 'Results\Original_data.xlsx';
P_data      = readtable(input_file,'Sheet','data');%,'Range','1:12',
tspan       = table2array(P_data(1,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'}));
Lys_ratio   = transpose(table2array(P_data(2,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'})));
P_ratio     = transpose(table2array(P_data(3:end,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'})));
EtaP        = transpose(table2array(P_data(3:end,'Lys_Conc_microM') ));

Lys_conc    = 206;%microM 
sum_EtaP    = 183190; %micoM
len_P_data  = (size(P_ratio,2)); 

%% Global optimization  using Firstst set of Ansatz.

Indi_fit_res_tab    = readtable('Results\Indi_fit_res.xlsx','Sheet','Indi_fit');
%Ansatz  = [ median(table2array(Indi_fit_res_tab(:,{'gamma_Lys'}))),transpose(table2array(Indi_fit_res_tab(:,{'gamma_Pi'}))), median(table2array(Indi_fit_res_tab(:,{'gamma_restP'})))];
deri_fit_res_tab    = readtable('Results\Derivative_fit_res.xlsx','Sheet','gamma_DeriFit');
Ansatz  = [ transpose(table2array(deri_fit_res_tab(1:end-1,{'gamma_deriFit'}))), median(table2array(Indi_fit_res_tab(:,{'gamma_restP'})))];
number_param    = 2 + len_P_data;
%FirstSetAnsatz = rand(1,number_param);
EtaP_temp       = [EtaP, sum_EtaP-sum(EtaP)];
lb              = zeros(1,number_param); 
ub              = []; 

%Global Optimization using Global search
opts        = optimoptions(@fmincon,'Algorithm','interior-point','UseParallel','always');%,
problem     = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)Glob_Opti_Lys_P(param, tspan,  Lys_ratio, P_ratio,  Lys_conc, EtaP_temp),'lb',lb,'ub',ub,'options',opts);
gs          = GlobalSearch('OutputFcn',@outputFunction_2,'Display','iter','NumTrialPoints',250);%
[Glob_fit_res,fval,exitflag,output,solutions] = run(gs,problem)

% Store the parameters 
Prot_name = [];
Prot_name{1}   = 'Lys';
for k = 1:1:len_P_data
    Prot_name{k+1} = sprintf('Protein_%d',k);
end
Prot_name{k+2}     = 'UIP';
Glob_fit_res_tab    = array2table(transpose(Glob_fit_res),'VariableNames',{'gamma'});
gamma_name_tab      = cell2table(transpose(Prot_name),'VariableNames',{'name'});
writetable([gamma_name_tab Glob_fit_res_tab],'Results\Glob_fit_res.xlsx','Sheet','Glob_fit');

exit
