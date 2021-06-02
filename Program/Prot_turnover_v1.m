function [] = Prot_turnover_v1(data,parameters)
%program to study the protein turnover kinetics. version 1.0.0
 
%read SILAC data 
size_Protein_data = size(data.SILAC_data,2);
fun_deg = @(b,tspan) exp(-b*tspan); 
for i =1:size_Protein_data
    Pdata_temp = data.SILAC_data(:,i);
	non_NaN_index_of_Pdata_temp = find(~isnan(Pdata_temp));
	Pdata_temp = Pdata_temp(non_NaN_index_of_Pdata_temp);
    rng default % For reproducibility
	ExpDegRate_temp = fminsearch(@(b) norm(Pdata_temp - fun_deg(b,data.t(non_NaN_index_of_Pdata_temp))), rand(1,1)) ;
	apparent_half_lives(i) = ExpDegRate_temp;
end
apparent_half_lives = transpose(log(2)./apparent_half_lives);
[~, idx] = sort(apparent_half_lives);
data.SILAC_data = data.SILAC_data(:,idx); 
data.ProtConc = data.ProtConc(:,idx);
data.protInfo = data.protInfo(idx, :); 

data.SILAC_data2 = data.SILAC_data(:, (sum(~isnan(data.SILAC_data)) > 4));
data.protInfo2 = data.protInfo((sum(~isnan(data.SILAC_data)) > 4),:);
data.ProtConc2 = data.ProtConc(:, (sum(~isnan(data.SILAC_data)) > 4));

% define Bin_size and Bin#
if parameters.Num_prot_to_fit < parameters.prot_per_fit
    Bin_size = parameters.Num_prot_to_fit;
else
    Bin_size = parameters.prot_per_fit;
end

Bins = ceil(size_Protein_data/Bin_size);
indices = ([0:Bin_size-1]*Bins)+1;
indices_0 = indices;
for i = 1:Bins-1
    indices = [indices (indices_0+i)];
end
indices(indices > size_Protein_data) = [];
data.SILAC_data = data.SILAC_data(:,indices); 
data.ProtConc = data.ProtConc(:,indices);
data.protInfo = data.protInfo(indices, :); 

data.SILAC_data = data.SILAC_data(:, 1:parameters.Num_prot_to_fit); 
data.ProtConc = data.ProtConc(:, 1:parameters.Num_prot_to_fit);
data.protInfo = data.protInfo(1:parameters.Num_prot_to_fit,:); 
size_Protein_data = size(data.SILAC_data,2);
%creating an array with bin size
se = [0];
for i = 1:floor(size_Protein_data/parameters.prot_per_fit)
    se(i+1) = i*parameters.prot_per_fit;
end
if (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) > 50)
    se(end+1) = size_Protein_data; 
elseif (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) < 51) && (length(se) >1)
    se(end) = size_Protein_data; 
elseif (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) < 51) && (length(se) == 1)
    se(end+1) = size_Protein_data; 
end

%% create araay and bin for selcted protein with all time points 
size_Protein_data2 = size(data.SILAC_data2,2);
% define Bin_size and Bin#
if parameters.Num_prot_to_fit < parameters.prot_per_fit
    Bin_size = parameters.Num_prot_to_fit;
else
    Bin_size = parameters.prot_per_fit;
end

Bins = ceil(size_Protein_data2/Bin_size);
indices = ([0:Bin_size-1]*Bins)+1;
indices_0 = indices;
for i = 1:Bins-1
    indices = [indices (indices_0+i)];
end
indices(indices > size_Protein_data2) = [];
data.SILAC_data2 = data.SILAC_data2(:,indices); 
data.ProtConc2 = data.ProtConc2(:,indices);
data.protInfo2 = data.protInfo2(indices, :);

if parameters.Num_prot_to_fit < size_Protein_data2
    data.SILAC_data2 = data.SILAC_data2(:, 1:parameters.Num_prot_to_fit); 
    data.ProtConc2 = data.ProtConc2(:, 1:parameters.Num_prot_to_fit);
    data.protInfo2 = data.protInfo2(1:parameters.Num_prot_to_fit, :); 
end
size_Protein_data2 = size(data.SILAC_data2,2);
%creating an array with bin size
se2 = [0];
for i = 1:floor(size_Protein_data2/parameters.prot_per_fit)
    se2(i+1) = i*parameters.prot_per_fit;
end
if (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) > 50)
    se2(end+1) = size_Protein_data2; 
elseif (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) < 51) && (length(se2) >1)
    se2(end) = size_Protein_data2; 
elseif (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) < 51) && (length(se2) == 1)
    se2(end+1) = size_Protein_data2; 
end

if length(se2) > 5
    se2 = se2(1:5);
end 

opts = optimset('FinDiffRelStep',1e-3);
opts2           = optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',1e-3,'FiniteDifferenceType','central', 'MaxIterations', 2, 'MaxFunctionEvaluations', 2);%,
ms              = MultiStart('UseParallel','always','Display','iter');
rng default % For reproducibility
Ansatz_          = rand(11,size_Protein_data+2);

%% Setting-3
if parameters.setting3 == 1
    % Fitting/optimizing data with all time points to get the Lys degredation rate 
    Error_free_Lys =  100; gama_Lys_best =  100;  prot_fitting = [];free_Lys_fitting = [];
    for i = 1:length(se2)-1
        fprintf('Optimizing proteins-%d to %d, out of %d  with setting-3\n',se2(i)+1,se2(i+1), se2(end))
        data.SILAC_data_temp = [data.Lys' data.SILAC_data2(:,se2(i)+1:se2(i+1))];
        data.Conc_temp  = [data.LysConc, data.ProtConc2(se2(i)+1:se2(i+1)), parameters.tot_Lys-sum(data.ProtConc2(se2(i)+1:se2(i+1)))];
        number_param    = 1 + size(data.SILAC_data_temp,2);%
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [1 (se2(i)+2:se2(i+1)+1) se2(end)+2]); 
        problem = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, 3),'lb',lb,'ub',ub,'options',opts);
        Ansatz          = Ansatz_(2:end, [1 (se2(i)+2:se2(i+1)+1) se2(end)+2]);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
        
        for j =1 : length(solutions)
            Lys_P_init          = ones(1,length(data.Conc_temp));
            [t,Lys_P]           = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,solutions(j).X,data.Conc_temp),data.t,Lys_P_init);
            Error_free_Lys_     = sum(abs(Lys_P(:,1) - data.Lys'));
            Error_proteins_Lys  = sum(nansum((data.SILAC_data_temp(:,1:end) - Lys_P(:,1:end-1)).^2));
            temp = [2; 2; 2; 2];
            if (length(temp(Lys_P(2:end,1) <= (data.Lys(2:end)'))) > 2) && (Error_free_Lys_ < Error_free_Lys ) %  -Error_free_Lys*0.1      && (solutions(j).fval < fval_best)
                Error_free_Lys  = Error_free_Lys_; 
             	gama_Lys_best   =  solutions(j).X(1);
            	Gfit            = solutions(j).X;
            end
        end
    end
    
    % Fitting/optimizing all data using the common paranters obtained in the above step
    Glob_fit_Prot = []; Glob_fit_Lys = [];  Glob_fit_UIP = []; res_error_ = [];
    for i = 1:length(se)-1
        fprintf('\n Fitting/Optimizing proteins-%d to %d, out of %d with setting-3 \n',se(i)+1,se(i+1), se(end))
        data.SILAC_data_temp = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
        number_param    = size(data.SILAC_data_temp,2);%
        data.Conc_temp  = [data.LysConc, data.ProtConc(se(i)+1:se(i+1)), parameters.tot_Lys-sum(data.ProtConc(se(i)+1:se(i+1)))];
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [(se(i)+2:se(i+1)+1) se(end)+2]);
        problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp,3, gama_Lys_best),'lb',lb,'ub',ub,'options',opts);
        Ansatz          = Ansatz_(2:end, [(se(i)+2:se(i+1)+1) se(end)+2]); 
        custpts         = CustomStartPointSet(Ansatz);
        [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
       
        Error_free_Lys = 100;
        for j =1 : length(solutions)
            Lys_P_init          = ones(1,length(data.Conc_temp));
            [t,Lys_P]           = ode15s(@(t,Lys_P)PT_ODE_fixLys(t,Lys_P,solutions(j).X, data.Conc_temp, gama_Lys_best),data.t,Lys_P_init);
            Error_free_Lys_     = sum(abs(Lys_P(:,1) - data.Lys'));
            Error_proteins_Lys  = sum(nansum((data.SILAC_data_temp(:,1:end) - Lys_P(:,1:end-1)).^2));
            if (length(temp(Lys_P(2:end,1) <= (data.Lys(2:end)'))) > 2) && (Error_free_Lys_ < Error_free_Lys)
                Error_free_Lys = Error_free_Lys_;
                Gfit = solutions(j).X;
            end
        end
        
        if Gfit(end) == 0; Gfit(end)  =1e-32; end;
        param_temp      = [log(2)./Gfit] ;
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2,param_temp,lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, 3, gama_Lys_best);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        
        %calculating the residual error
        Lys_P_init      = ones(1,length(data.Conc_temp));
        [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_fixLys(t,Lys_P,log(2)./Opts, data.Conc_temp, gama_Lys_best),data.t, Lys_P_init);
        res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,1:end-1)).^2));
        res_error_(i,:) = [sum(res_error(2:end)) nansum((data.Lys' - Lys_P(:,1)).^2)];
        fprintf('\n residual error  %6.4f \n', sum(res_error))
        Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(1:end-1)',(((Ci(1:end-1,2)-Ci(1:end-1,1))/2)/2),res_error(2:end)'];%[param_temp(2:end-1)', res_error(2:end)'];
        Glob_fit_Lys(i,:) = log(2)/gama_Lys_best;
        Glob_fit_UIP(i,:) = param_temp(end);
    end
    fprintf('\n *****  Completed Optimization of paramters ******** \n \n Now exporting half-live to excel file \n')
    % Store the parameters 
    synthesis_rate = (data.ProtConc)'.* (log(2)./Glob_fit_Prot(:,1));
    Gfit_tab    = array2table([ data.SILAC_data', Glob_fit_Prot,  synthesis_rate, data.ProtConc' ],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI','residual_error', 'synthesis_rate', 'protein_conc'});
    writetable([data.protInfo, Gfit_tab],parameters.out_file, 'Sheet','Setting3');
end

%% Setting-2
if parameters.setting2 == 1
    % Fitting/optimizing data with all time points to get the Lys degredation rate and ratio between avearage protein-bound Lys and free lys
    Error_free_Lys =  100; gama_Lys_best =  100; k=0; free_LysError = 0; prot_fitting = [];free_Lys_fitting = [];
    while (gama_Lys_best == 100) && (free_LysError <= 50)
        for i = 1:length(se2)-1
            fprintf('Optimizing proteins-%d to %d, out of %d  with setting-2\n',se2(i)+1,se2(i+1), se2(end))
            data.SILAC_data_temp = [(data.Lys'-data.Lys'*free_LysError/100) data.SILAC_data2(:,se2(i)+1:se2(i+1))];
            number_param    = 1 + size(data.SILAC_data_temp,2);%
            lb              = zeros(1,number_param); 
            ub              = [];
            Ansatz          = Ansatz_(1, [1 (se2(i)+2:se2(i+1)+1) se2(end)+2]); 
            problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,0,2),'lb',lb,'ub',ub,'options',opts);
            Ansatz          = Ansatz_(2:end, [1 (se2(i)+2:se2(i+1)+1) se2(end)+2]);
            custpts = CustomStartPointSet(Ansatz);
            [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);

            for j =1 : length(solutions)
                Lys_P_init      = ones(1,size(data.SILAC_data_temp,2));
                [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t,Lys_P_init);
                Error_free_Lys_ = sum(abs(Lys_P(:,1) - (data.Lys')));
                Error_proteins_Lys = sum(nansum((data.SILAC_data_temp(:,1:end) - Lys_P(:,1:end)).^2));
                temp = [2; 2; 2; 2];
                if (length(temp(Lys_P(2:end,1) <= data.Lys(2:end)')) > 2) && (Error_free_Lys_ < Error_free_Lys) 
                    Error_free_Lys = Error_free_Lys_; 
                    gama_Lys_best =  solutions(j).X(1);
                    Gfit = solutions(j).X;
                end
            end
        end
        k = k+1; free_LysError = 10*k;
    end
    
    % Fitting/optimizing all data using the common paranters obtained in the above step
    Glob_fit_Prot = []; Glob_fit_Lys = []; Glob_fit_ratioConst = []; res_error_ = []; 
    for i = 1:length(se)-1
        fprintf('\nOptimizing proteins-%d to %d, out of %d  with setting-2\n',se(i)+1,se(i+1), se(end))
        data.SILAC_data_temp = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
        number_param    = size(data.SILAC_data_temp,2);%
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [(se(i)+2:se(i+1)+1), se2(end)+2]); 
        problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,0,2, gama_Lys_best),'lb',lb,'ub',ub,'options',opts);
        Ansatz          = Ansatz_(2:end, [(se(i)+2:se(i+1)+1), se2(end)+2]);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
        
        Error_free_Lys = 100;
        for j =1 :length(solutions)
            Lys_P_init      = ones(1,size(data.SILAC_data_temp,2));
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,solutions(j).X, gama_Lys_best),data.t,Lys_P_init);
            Error_free_Lys_ = sum(abs(Lys_P(2:end,1) - data.Lys(2:end)'));
            temp = [2; 2; 2; 2];
            if (length(temp(Lys_P(2:end,1) <= (data.Lys(2:end)'))) > 2) && (Error_free_Lys_ < Error_free_Lys )
                Error_free_Lys = Error_free_Lys_;
                Gfit = solutions(j).X;
            end
        end
              
        param_temp      = [log(2)./Gfit(1:end-1), Gfit(end)] ;
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2, param_temp, lb, ub,opts2, data.t, data.SILAC_data_temp, 0, 2,gama_Lys_best);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        %calculating the residual error
        Lys_P_init      = ones(1,size(data.SILAC_data_temp,2));
        [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts(1:end-1), Opts(end)], gama_Lys_best),data.t,Lys_P_init);
        res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,1:end)).^2));
        res_error_(i,:) = [sum(res_error(2:end)) nansum((data.Lys' - Lys_P(:,1)).^2)];
        fprintf( " LSE for the solustion choosen  %4f ", sum(res_error_(i,:)));
        Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(1:end-1)', ((Ci(1:end-1,2)-Ci(1:end-1,1))/2)/2,res_error(2:end)'];
    end
    fprintf('\n *****  Completed Optimization of paramters ******** \n Now exporting half-live to excel file \n')
    Gfit_tab    = array2table([ data.SILAC_data', Glob_fit_Prot],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI','residual_error'});
    writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting2');
end

%% Setting-1
if parameters.setting1 == 1
    % Fitting/optimizing data with all time points to get the Lys degredation rate and ratio between avearage protein-bound Lys and free lys
    Error_free_Lys =  100;     gama_Lys_best =  100; k = 1; free_LysError = 10;
    while (gama_Lys_best == 100) && (free_LysError <= 50)
        for i = 1:length(se2)-1
            fprintf('Optimizing proteins-%d to %d, out of %d  with setting-1\n',se2(i)+1,se2(i+1), se2(end))
            data.SILAC_data_temp = data.SILAC_data2(:,se2(i)+1:se2(i+1));
            P_ratio_fastProt = [min(data.SILAC_data2(2,:)) min(data.SILAC_data2(3,:)) min(data.SILAC_data2(4,:)) min(data.SILAC_data2(5,:))];
            data.SILAC_data_temp = [[1 (P_ratio_fastProt - P_ratio_fastProt*(free_LysError/100))]' data.SILAC_data_temp];
            number_param    =1 + size(data.SILAC_data_temp,2);%
            lb              = zeros(1,number_param); 
            ub              = [];
            Ansatz          = Ansatz_(1, [1 (se2(i)+2:se2(i+1)+1) se2(end)+2]); 
            problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,0,2),'lb',lb,'ub',ub,'options',opts);
            ms              = MultiStart('UseParallel','always','Display','iter');
            Ansatz          = Ansatz_(2:end, [1 (se2(i)+2:se2(i+1)+1) se2(end)+2]);
            custpts = CustomStartPointSet(Ansatz);
            [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);

            for j =1 : length(solutions)
                Lys_P_init      = ones(1, size(data.SILAC_data_temp,2));%Lys_P_init      = ones(1,1 + size(data.SILAC_data_temp,2));
                [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t,Lys_P_init);
                Error_free_Lys_ = sum(abs(Lys_P(2:end,1) - P_ratio_fastProt'));
                temp = [2; 2; 2; 2];
                if (length(temp(Lys_P(2:end,1) < P_ratio_fastProt')) > 2) && (Error_free_Lys_ < Error_free_Lys)
                    Error_free_Lys  = Error_free_Lys_;
                    Gfit            = solutions(j).X;
                    gama_Lys_best   =  solutions(j).X(1);
                end 
            end
        end
        k= k+1; free_LysError = 10*k;
    end
    
    % Fitting/optimizing all data using the common paranters obtained in the above step
    Glob_fit_Prot = []; Glob_fit_Lys = []; Glob_fit_ratioConst = []; res_error_ = []; res_error_Lys = [];
    for i = 1:length(se)-1
        fprintf('Optimizing proteins-%d to %d, out of %d with setting-1\n',se(i)+1,se(i+1), se(end))
        data.SILAC_data_temp = data.SILAC_data(:,se(i)+1:se(i+1));
        number_param    = 1+ size(data.SILAC_data_temp,2);%
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [(se(i)+2:se(i+1)+1), se(end)+2]); 
        problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,0,1,gama_Lys_best ),'lb',lb,'ub',ub,'options',opts);
        ms              = MultiStart('UseParallel','always','Display','iter');
        Ansatz          = Ansatz_(2:end, [(se(i)+2:se(i+1)+1), se(end)+2]);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,~,~,~,~]  = run(ms,problem,custpts); 
        
        param_temp      = [log(2)./Gfit(1:end-1), Gfit(end)] ;
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2,param_temp,lb,ub,opts2, data.t, data.SILAC_data_temp,  0, 1,gama_Lys_best);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        Lys_P_init      = ones(1,1 + size(data.SILAC_data_temp,2));
        [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts(1:end)], gama_Lys_best),data.t,Lys_P_init);
        res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,2:end)).^2));
        Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(1:end-1)',((Ci(1:end-1,2)-Ci(1:end-1,1))/2)/2,res_error(1:end)'];%[param_temp(2:end-1)',res_error(1:end)'];
        Glob_fit_Lys(i,:) = log(2)/gama_Lys_best;
        Glob_fit_ratioConst(i,:) = Opts(end);
    end
    fprintf('\n *****  Completed Optimization of paramters with setting-1 ********\n Now exporting half-live to excel file \n')
    Gfit_tab    = array2table([ data.SILAC_data', Glob_fit_Prot ],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI','residual_error'});
    writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting1');
end
fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')
end

