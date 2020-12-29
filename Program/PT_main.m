%% program to study the protein turnover kinetics. 5/20/2019

close all
clearvars 

%% read the original data
input_file  = 'test_data.xlsx';
out_file = 'Half-lives.xlsx';
sum_Conc    = 41300; %micoM
d.t_long = linspace(0, 32,321);

All_data   = readtable(input_file);
data.t     = table2array(All_data(1,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'}));  
data.Lys   = table2array(All_data(2,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'}));% pulse times 
data.Conc  = table2array(All_data(2:end,{'Concentration_microM'}))';% pulse times 
data.t_long = linspace(0, 32,321);
data.protInfo  = table2array(All_data(3:end,{'Name'}));                       % protein IDs

SILAC_data      = table2cell(All_data(3:end,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'}));          % data range
SILAC_data_NaN  = cellfun(@(x) ~isa(x,'double'),SILAC_data);     % check for possible 0x0 char reads
SILAC_data(SILAC_data_NaN)               = {NaN};                                % set those to NaNs
data.SILAC_data                 = cell2mat(SILAC_data)';                          % heavy/light ratios,



Num_prot_for_optim = 100;
se = [0];
for i = 1:floor(length(data.SILAC_data)/Num_prot_for_optim)
    se(i+1) = i*Num_prot_for_optim;
end
if mod(length(data.SILAC_data),Num_prot_for_optim) > 0; se(end+1) = length(data.SILAC_data); end

Glob_fit_Prot = []; Glob_fit_Lys = [];  Glob_fit_UIP = [];
for i = 1:length(se)-1
    fprintf('Optimizing proteins-%d to %d, out of %d \n',se(i)+1,se(i+1), se(end))
    data.SILAC_data_temp = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
    number_param    = 1 + size(data.SILAC_data_temp,2);%
    data.Conc_temp     = [data.Conc(1), data.Conc(se(i)+2:se(i+1)+1), sum_Conc-sum(data.Conc(se(i)+2:se(i+1)+1))];
    lb              = zeros(1,number_param); 
    ub              = [];
    Ansatz          = rand(1,number_param);
    opts            = optimoptions(@fmincon,'Algorithm','sqp','UseParallel','always');%,
    problem         = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO(param, data.t,  data.SILAC_data_temp,  data.Conc_temp),'lb',lb,'ub',ub,'options',opts);

    ms              = MultiStart('UseParallel','always','Display','iter');
    Ansatz          = rand(10,number_param); custpts = CustomStartPointSet(Ansatz);
    [Gfit,~,~,~,~]  = run(ms,problem,custpts);    
    
    param_temp      = (log(2)./Gfit);if param_temp(end) == Inf; param_temp(end)  =log(2)/1e-16; end; 
    [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin,param_temp,lb,ub,[], data.t, data.SILAC_data_temp,  data.Conc_temp);
    Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
    %calculating the residual error
    Lys_P_init      = ones(1,length(data.Conc_temp));
    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,log(2)./Opts, data.Conc_temp),data.t,Lys_P_init);
    res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,1:end-1)).^2));
    Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(2:end-1)',Ci(2:end-1,:),res_error(2:end)'];
end

fprintf('\n *****  Completed Optimization of paramters; Now exporting half-live to excel file **** \n')

% Store the parameters
N  = size(All_data,1)-1;
data_table = cell([N 5]);
data_table(1, 1:5) = {'Protein_Description', 'HalfLife(days)', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error'};
data_table(2:end,1) = data.protInfo;
    
%Tc(1,1+(1:N)) = data.cnd(data.cndI)
data_table(2:end,2:end) = num2cell(Glob_fit_Prot); %num2cell(data.t)
data_table = cell2table(data_table);
writetable(data_table,out_file,'WriteVariableNames',false);% 'Sheet','half-lives'

fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')

