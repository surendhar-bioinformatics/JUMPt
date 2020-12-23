%% program to study the protein turnover kinetics. 5/20/2019

close all
clearvars 

%% read the original data
input_file  = 'test_data.xlsx';
[~,~,raw]   = xlsread(input_file);

cnd                 = raw(1,2:end);
cnd = cnd(cellfun('isclass', cnd, 'char'));
N                   = numel(cnd)+1;
[d.cnd, ~, d.cndI]  = unique(cnd,'stable'); d.cndI=d.cndI'; % unique conditions and corresponding data indices
d.t                 = cell2mat(raw(2,4:N-1));  
d.Lys                 = cell2mat(raw(3,4:N-1));% pulse times        
d.t_long = linspace(0, 32,321);
d.Conc              = cell2mat(raw(3:end,N))';                 % pulse times        
d.pID               = raw(4:end,1);                         % protein IDs
Y                   = raw(3+(1:numel(d.pID)),4:N-1);          % data range
Yc                  = cellfun(@(x) ~isa(x,'double'),Y);     % check for possible 0x0 char reads
Y(Yc)               = {NaN};                                % set those to NaNs
d.y                 = cell2mat(Y)';                          % heavy/light ratios,
sum_Conc    = 41300; %micoM
d.t_long = linspace(0, 32,321);
out_file = 'Half-lives.xlsx';

%% Global optimization  using Firstst set of Ansatz.

%% Global Optimization using Multi search
Num_prot_for_optim = 100;
se = [0];
for i = 1:floor(length(d.y)/Num_prot_for_optim)
    se(i+1) = i*Num_prot_for_optim;
end
if mod(length(d.y),Num_prot_for_optim) > 0; se(end+1) = length(d.y); end

Glob_fit_Prot = []; Glob_fit_Lys = [];  Glob_fit_UIP = [];
for i = 1:length(se)-1
    fprintf('Optimizing proteins-%d to %d, out of %d \n',se(i)+1,se(i+1), se(end))
    d.y_temp = [d.Lys' d.y(:,se(i)+1:se(i+1))];
    number_param    = 1 + size(d.y_temp,2);%
    d.Conc_temp     = [d.Conc(1), d.Conc(se(i)+2:se(i+1)+1), sum_Conc-sum(d.Conc(se(i)+2:se(i+1)+1))];%d.Conc_temp       = [d.Conc];
    lb              = zeros(1,number_param);%+0.009 
    ub              = [];
    Ansatz          = rand(1,number_param);
    opts            = optimoptions(@fmincon,'Algorithm','sqp','UseParallel','always', 'MaxIterations', 50000, 'MaxFunctionEvaluations', 500000000);%,
    problem         = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO(param, d.t,  d.y_temp,  d.Conc_temp),'lb',lb,'ub',ub,'options',opts);

    ms              = MultiStart('UseParallel','always','Display','iter');
    Ansatz          = rand(10,number_param); custpts = CustomStartPointSet(Ansatz);
    [Gfit,~,~,~,~]  = run(ms,problem,custpts);    
    
    param_temp      = (log(2)./Gfit);if param_temp(end) == Inf; param_temp(end)  =log(2)/1e-16; end; 
    [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin,param_temp,lb,ub,[], d.t, d.y_temp,  d.Conc_temp);
    Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
    %calculating the residual error
    Lys_P_init      = ones(1,length(d.Conc_temp));
    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,log(2)./Opts, d.Conc_temp),d.t,Lys_P_init);
    res_error       = (nansum((d.y_temp - Lys_P(:,1:end-1)).^2));
    Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(2:end-1)',Ci(2:end-1,:),res_error(2:end)'];
end

fprintf('\n *****  Completed Optimization of paramters; Now exporting half-live to excel file **** \n')

% Store the parameters 
Gfit_tab    = array2table(Glob_fit_Prot,'VariableNames',{ 'HalfLife', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error' });
writetable(Gfit_tab,out_file, 'Sheet','half-lives');


fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')


