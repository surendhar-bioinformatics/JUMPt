function [] = Prot_turnover(data,parameters)
%program to study the protein turnover kinetics. 5/20/2019
%   Detailed explanation goes here
% 
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
% define Bin_size and Bin#
if size_Protein_data < parameters.prot_per_fit
    Bin_size = size_Protein_data;
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

%creating an array with bin size
se = [0];
for i = 1:floor(size(data.SILAC_data,2)/parameters.prot_per_fit)
    se(i+1) = i*parameters.prot_per_fit;
end
if (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 0)  && (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 50)
    se(end+1) = size(data.SILAC_data,2); 
elseif (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 0)  && (mod(size(data.SILAC_data,2),parameters.prot_per_fit) < 51) && (length(se) >1)
    se(end) = size(data.SILAC_data,2); 
elseif (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 0)  && (mod(size(data.SILAC_data,2),parameters.prot_per_fit) < 51) && (length(se) == 1)
    se(end+1) = size(data.SILAC_data,2); 
end

opts            = optimoptions(@fmincon,'Algorithm','sqp','UseParallel','always','MaxIterations', 50000000, 'MaxFunctionEvaluations', 50000000000);%,
rng default % For reproducibility
Ansatz_          = rand(11,size_Protein_data+2);

%Binning all the proteins
se = [0];
for i = 1:floor(size(data.SILAC_data,2)/parameters.prot_per_fit)
    se(i+1) = i*parameters.prot_per_fit;
end
if (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 0)  && (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 50)
    se(end+1) = size(data.SILAC_data,2); 
elseif (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 0)  && (mod(size(data.SILAC_data,2),parameters.prot_per_fit) < 51) && (length(se) >1)
    se(end) = size(data.SILAC_data,2); 
elseif (mod(size(data.SILAC_data,2),parameters.prot_per_fit) > 0)  && (mod(size(data.SILAC_data,2),parameters.prot_per_fit) < 51) && (length(se) == 1)
    se(end+1) = size(data.SILAC_data,2); 
end

%% Setting-3
if parameters.setting3 == 1
    Glob_fit_Prot = []; Glob_fit_Lys = [];  Glob_fit_UIP = []; 
    for i = 1:length(se)-1
        fprintf('Optimizing proteins-%d to %d, out of %d \n',se(i)+1,se(i+1), se(end))
        data.SILAC_data_temp = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
        number_param    = 1 + size(data.SILAC_data_temp,2);%
        data.Conc_temp  = [data.LysConc, data.ProtConc(se(i)+1:se(i+1)), parameters.tot_Lys-sum(data.ProtConc(se(i)+1:se(i+1)))];
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [1 (se(i)+2:se(i+1)+1) se(end)+2]);
        problem = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO(param, data.t,  data.SILAC_data_temp,data.Conc_temp,3),'lb',lb,'ub',ub,'options',opts);
        ms              = MultiStart('UseParallel','always','Display','iter');
        Ansatz          = Ansatz_(2:end, [1 (se(i)+2:se(i+1)+1) se(end)+2]); 
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,fval,exitflag,output,solutions]  = run(ms,problem,custpts);
        P_ratio_fastProt = data.Lys(2:end);
        Error_free_Lys = 1000;
        for j =1 : length(solutions)
            Lys_P_init      = ones(1,length(data.Conc_temp));
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,solutions(j).X, data.Conc_temp),data.t,Lys_P_init);
            Error_free_Lys_ = sum(abs(Lys_P(2:end,1) - P_ratio_fastProt'));
            if Error_free_Lys_ < Error_free_Lys
                Error_free_Lys = Error_free_Lys_;
                Gfit = solutions(j).X;
            end
        end
        if Gfit(end) == 0; Gfit(end)  =1e-32; end;
        param_temp      = [log(2)./Gfit] ;%if param_temp(end) == Inf; param_temp(end)  =log(2)/1e-16; end;
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin,param_temp,lb,ub,[], data.t, data.SILAC_data_temp,  data.Conc_temp, 3);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        %calculating the residual error
        Lys_P_init      = ones(1,length(data.Conc_temp));
        [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,log(2)./Opts, data.Conc_temp),data.t,Lys_P_init);
        res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,1:end-1)).^2));
        Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(2:end-1)',Ci(2:end-1,:),res_error(2:end)'];%[param_temp(2:end-1)', res_error(2:end)'];
        Glob_fit_Lys(i,:) = [param_temp(1)];
        Glob_fit_UIP(i,:) = param_temp(end);
     end
    Glob_fit_res0 = [median(Glob_fit_Lys(:,1)),Glob_fit_Prot(:,1)',median(Glob_fit_UIP(:,1))];
    Glob_fit_res = log(2)./Glob_fit_res0;
    data.SILAC_data_temp = [data.Lys' data.SILAC_data];
    data.Conc_temp  = [data.LysConc, data.ProtConc, parameters.tot_Lys-sum(data.ProtConc)];
    fprintf('\n *****  Completed Optimization of paramters ******** \n \n Now exporting half-live to excel file \n')
    % Store the parameters 
    Gfit_tab    = array2table(Glob_fit_Prot,'VariableNames',{'HalfLife_in_days', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error'});
    writetable([data.protInfo, Gfit_tab],parameters.out_file, 'Sheet','Setting3');
    
    % plot the Raw data
    figure(); hold on; set(gca,'LineWidth',2.5); set(gca,'FontSize',16); % box on; grid on;
    plot(data.t, data.Lys,'s-r','MarkerSize',9,'LineWidth',2);
    plot(data.t, data.SILAC_data,'s-b','MarkerSize',8,'LineWidth',2);
    xlabel('Time (days)', 'FontSize',16); ylabel('Fraction of  ''light''  Lys ','FontSize',16);
    title('Raw data', 'FontSize',16);
    hlegend= legend('free-Lys','Proteins'); set(hlegend, 'FontSize',16, 'Box','off'); ylim([0 1.03]); xlim([0 33]); hold off
    saveas(gcf,'Raw-data.png')

    % plot the Fitting results
    Lys_P_init      = ones(1,length(data.Conc_temp));
    [t,Lys_P_simu]  = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,Glob_fit_res, data.Conc_temp),data.t_long,Lys_P_init);
    figure();hold on; set(gca,'LineWidth',2.5); set(gca,'FontSize',16); % box on; grid on;
    plot(t, Lys_P_simu(:,1),'-r','MarkerSize',8,'MarkerEdgeColor','g','LineWidth',2.5);
    plot(t, Lys_P_simu(:,2),'-b','MarkerSize',8,'LineWidth',2);
    plot(t, Lys_P_simu(:,end),'-c','MarkerSize',8,'LineWidth',2);
    plot(t, Lys_P_simu(:,2:end-1),'-b','MarkerSize',8,'LineWidth',2);
    plot(data.t, data.SILAC_data,'sb','MarkerSize',8,'LineWidth',2);
    plot(data.t, data.Lys,'sr','MarkerSize',9,'LineWidth',2);
    plot(t, Lys_P_simu(:,1),'--r','MarkerSize',8,'MarkerEdgeColor','g','LineWidth',2.5);
    xlabel('Time (days)', 'FontSize',16); ylabel('Fraction of  ''light''  Lys ','FontSize',16);
    title('proteins and Lys fitting(Setting-3)', 'FontSize',16);
    hlegend= legend('free-Lys','Proteins','UIP'); set(hlegend, 'FontSize',16, 'Box','off'); ylim([0 1.03]); xlim([0 33]); hold off
    saveas(gcf,'Setting-3.png')
    %fprintf('\n *****  Completed exporting half-live to excel file ******** \n ******Program is complete*******\n\n')
end

%% Setting-2
if parameters.setting2 == 1
    Glob_fit_Prot = []; 
    for i = 1:length(se)-1
        fprintf('Optimizing proteins-%d to %d, out of %d  with setting-2\n',se(i)+1,se(i+1), se(end))
        data.SILAC_data_temp = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
        number_param    = 1 + size(data.SILAC_data_temp,2);%
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [1 (se(i)+2:se(i+1)+1) se(end)+2]); 
        problem         = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO(param, data.t, data.SILAC_data_temp, 0,2),'lb',lb,'ub',ub,'options',opts);
        ms              = MultiStart('UseParallel','always','Display','iter');
        Ansatz          = Ansatz_(2:end, [1 (se(i)+2:se(i+1)+1) se(end)+2]);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,fval,exitflag,output,solutions]  = run(ms,problem,custpts);
        Gfit_orig = Gfit;
        P_ratio_fastProt = data.Lys(2:end);
        Error_free_Lys = 100;
        for j =1 : length(solutions)
            Lys_P_init      = ones(1,size(data.SILAC_data_temp,2));
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t,Lys_P_init);
            Error_free_Lys_ = sum(abs(Lys_P(2:end,1) - P_ratio_fastProt'));
            if Error_free_Lys_ < Error_free_Lys
                Error_free_Lys = Error_free_Lys_;
                Gfit = solutions(j).X;
            end
        end
        param_temp      = [log(2)./Gfit(1:end-1),Gfit(end)] ;%if param_temp(end) == Inf; param_temp(end)  =log(2)/1e-16; end; 
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin,param_temp,lb,ub,[], data.t, data.SILAC_data_temp,  0, 2);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        %calculating the residual error
        Lys_P_init      = ones(1,size(data.SILAC_data_temp,2));
        [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,[log(2)./Opts(1:end-1), Opts(end)]),data.t,Lys_P_init);
        res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,1:end)).^2));
        Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(2:end-1)',Ci(2:end-1,:),res_error(2:end)'];%[param_temp(2:end-1)',res_error(2:end)'];
    end
    fprintf('\n *****  Completed Optimization of paramters ******** \n Now exporting half-live to excel file \n')
    Gfit_tab    = array2table(Glob_fit_Prot,'VariableNames',{'HalfLife_in_days', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error' });
    writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting2');
end

%% Setting-1
if parameters.setting1 == 1
    Glob_fit_Prot = []; 
    for i = 1:length(se)-1
        fprintf('Optimizing proteins-%d to %d, out of %d with setting-1\n',se(i)+1,se(i+1), se(end))
        data.SILAC_data_temp = data.SILAC_data(:,se(i)+1:se(i+1));
        number_param    = 2 + size(data.SILAC_data_temp,2);%
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = Ansatz_(1, [1 (se(i)+2:se(i+1)+1) se(end)+2]); 
        problem        = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO(param, data.t, data.SILAC_data_temp, 0,1),'lb',lb,'ub',ub,'options',opts);
        ms              = MultiStart('UseParallel','always','Display','iter');
        Ansatz          = Ansatz_(2:end, [1 (se(i)+2:se(i+1)+1) se(end)+2]);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,fval,exitflag,output,solutions]  = run(ms,problem,custpts); 
        param_temp      = [log(2)./Gfit(1:end-1),Gfit(end)] ;%if param_temp(end) == Inf; param_temp(end)  =log(2)/1e-16; end; 
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin,param_temp,lb,ub,[], data.t, data.SILAC_data_temp,  0, 1);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);%calculating the residual error
        Lys_P_init      = ones(1,1 + size(data.SILAC_data_temp,2));
        [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,[log(2)./Opts(1:end-1), Opts(end)]),data.t,Lys_P_init);
        res_error       = (nansum((data.SILAC_data_temp - Lys_P(:,2:end)).^2));
        Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(2:end-1)',Ci(2:end-1,:),res_error(1:end)'];%[param_temp(2:end-1)',res_error(1:end)'];
    end
    fprintf('\n *****  Completed Optimization of paramters with setting-1 ********\n Now exporting half-live to excel file \n')
    Gfit_tab    = array2table(Glob_fit_Prot,'VariableNames',{'HalfLife_in_days', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error' });
    writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting1');
end
fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')
end

