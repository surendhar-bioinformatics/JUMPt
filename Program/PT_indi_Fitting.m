%% program to study the protein turnover kinetics.
% Performing individual fitting of each protein along with Lys and UIP and
% calculates its derivatives

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
tspan_long  = linspace(0, 32,3201);
tspan_medium = linspace(0, 32,321);

%% Find the parameter estimates for proteins one by one  including UIP.
Indi_fit_res    = [];
f_Lys_Protein   = [];
df_Lys_Protein  = [];
for i = 1:len_P_data
    Pi_ratio        = P_ratio(:,i);
    number_param    = 3; % gamma_Lys, gamma_Pi, gamma_UIP
    param_guess 	= rand(1,number_param);
    lb              = zeros(1,number_param); 
    ub              = [];
    EtaP_temp       = [EtaP(i), sum_EtaP-EtaP(i)];

    %Global Optimization using Global search
    opts        = optimoptions(@fmincon,'Algorithm','interior-point','UseParallel','always');
    problem     = createOptimProblem('fmincon','x0',param_guess,'objective',@(param)Glob_Opti_P(param, tspan, Lys_ratio, Pi_ratio, Lys_conc, EtaP_temp),'lb',lb,'ub',ub,'options',opts);
    GlobSearch  = GlobalSearch('Display','iter','NumTrialPoints',200);%
    [param_indFit,fval,exitflag_noUIP,output,solutions] = run(GlobSearch,problem);
    
    %Storing the parameters to an array
    Indi_fit_res(i,1:3) = param_indFit;
    Indi_fit_res(i,4:6) = log(2)./param_indFit;
    %Simulate the data using individual optimized parameters
    Lys_P_init          = ones(1,3);
    [t,Lys_P_simu]      = ode15s(@(t,Lys_P)PT_ODE(t, Lys_P, param_indFit, Lys_conc, EtaP_temp),tspan, Lys_P_init);
    Indi_fit_res(i, 7)  = sum(nansum((Lys_ratio- Lys_P_simu(:,1)).^2)); %LSE_Lys
    Indi_fit_res(i, 8)  = sum(nansum((Pi_ratio- Lys_P_simu(:,2)).^2)); %LSE_Pi
    
    %calculate derivatives of Proteins
    [t,Lys_P_simu]      = ode15s(@(t,Lys_P)PT_ODE(t, Lys_P, param_indFit, Lys_conc, EtaP_temp),tspan_medium, Lys_P_init);
    f_Lys_Protein(i+1,:) = Lys_P_simu(:,2);
    df_Lys_Protein(i+1,:) =  gradient(Lys_P_simu(:,2))./gradient(t);
        
end %for i = 1:len_p_data
%calculate derivatives of Lys data
gamma_Lys           = Indi_fit_res(:,1);
%gamma_Lys_median_index  = find(gamma_Lys==median(gamma_Lys));
[temp gamma_Lys_median_index] = min(abs(gamma_Lys-median(gamma_Lys)));
[t,Lys_P_simu]      = ode15s(@(t,Lys_P)PT_ODE(t, Lys_P, Indi_fit_res(gamma_Lys_median_index,1:3), Lys_conc, EtaP_temp),tspan_medium, Lys_P_init);
f_Lys_Protein(1,:)  = Lys_P_simu(:,1);
df_Lys_Protein(1,:) =  gradient(Lys_P_simu(:,1))./gradient(t);

Indi_fit_res_tab = array2table(Indi_fit_res,'VariableNames',{'gamma_Lys', 'gamma_Pi', 'gamma_restP','halfLife_Lys', 'halfLife_Pi', 'halfLife_restP', 'LSE_Lys', 'LSE_Pi'});
writetable(Indi_fit_res_tab,'Results\Indi_fit_res.xlsx','Sheet','Indi_fit');

Prot_name = [];
Prot_name{1}   = 'Lys';
for k = 1:1:len_P_data
    Prot_name{k+1} = sprintf('Protein_%d',k);
end
f_Lys_Protein_tab = array2table(transpose(f_Lys_Protein),'VariableNames',Prot_name);
writetable(f_Lys_Protein_tab,'Results\Indi_fit_res.xlsx','Sheet','f_Lys_P');
df_Lys_Protein_tab = array2table(transpose(df_Lys_Protein),'VariableNames',Prot_name);
writetable(df_Lys_Protein_tab,'Results\Indi_fit_res.xlsx','Sheet','df_Lys_P');
exit
