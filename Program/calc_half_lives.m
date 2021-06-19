function [] = calc_half_lives(data,params)
% Calculates the the half-lives for based on settings
% developed by Surendhar Reddy Chepyala and Junmin Peng
% version (2.0.0) 
opts	= optimoptions(@fmincon,'Algorithm','sqp','UseParallel','always','MaxIterations', 50000, 'MaxFunctionEvaluations', 500000);%,
opts2	= optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',1e-3,'FiniteDifferenceType','central', 'MaxIterations', 2, 'MaxFunctionEvaluations', 2,'Display','off');%,
opts3   = optimset('FinDiffRelStep',1e-3);
ms      = MultiStart('UseParallel','always','Display','off');
outfile	= join(['results', params.input_file], "_");
SILAC_food_impurity = (100-params.purity)/100;

if params.setting == 1
 	SimLys_gteq_expLys = 0; free_LysError = 0; k = 0; gama_Lys_best =  100;  min_Error_proteins_only = 100; min_Error_prot= 100;  
    while (free_LysError <= 20 ) %  
        fprintf('\nOptimizing prot-%d to %d to find free-Lys data with setting-%d \n',data.protBins_allTimes(1)+1,data.protBins_allTimes(2), params.setting)
        number_param	= 2 + size(data.SILAC_data_allTimes(:,data.protBins_allTimes(1)+1:data.protBins_allTimes(2)),2);%
        lb              = zeros(1,number_param); ub = [];
        Ansatz          = data.Ansatz_(1, 1:number_param); 

    	data.SILAC_data_temp    = [data.SILAC_data_allTimes(:,data.protBins_allTimes(1)+1:data.protBins_allTimes(2))];
      	fast_data           = [min(data.SILAC_data_allTimes(2,:)) min(data.SILAC_data_allTimes(3,:)) min(data.SILAC_data_allTimes(4,:)) min(data.SILAC_data_allTimes(5,:))]';
    	data.SILAC_data_temp    = [[1; fast_data - fast_data*(free_LysError/100)] data.SILAC_data_temp];
     	data.Conc_temp      = 0;
        if params.opti_algo == 2
            problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, SILAC_food_impurity, params.setting+1),'lb',lb,'ub',ub,'options',opts3);
        else
        	problem      	 = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO_fmincon(param, data.t,  data.SILAC_data_temp,data.Conc_temp, SILAC_food_impurity, params.setting+1),'lb',lb,'ub',ub,'options',opts);
        end
        
        Ansatz	= data.Ansatz_(2:end, 1:number_param);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
        for j =1 : length(solutions)
            Lys_P_init      = ones(1,number_param-1);
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X, SILAC_food_impurity),data.t,Lys_P_init);
            gama_Lys_best   =  solutions(j).X(1);
         	Error_proteins_only = sum(nansum((data.SILAC_data_temp(:,2:end) - Lys_P(:,2:end)).^2));

            temp = [2; 2; 2; 2];
            if (length(temp(Lys_P(2:end,1) <= fast_data)) > SimLys_gteq_expLys) || ((length(temp(Lys_P(2:end,1) <= fast_data)) == SimLys_gteq_expLys ) && (Error_proteins_only < min_Error_prot )) 
                SimLys_gteq_expLys = length(temp(Lys_P(2:end,1) <= fast_data));
                PredLys         = Lys_P(:,1);
                min_Error_prot  = Error_proteins_only;
                gama_Lys_best   =  solutions(j).X(1);
                Concen_ratio    = solutions(j).X(end);
                Gfit_final      = solutions(j).X;
            elseif Error_proteins_only < min_Error_proteins_only
                min_Error_proteins_only = Error_proteins_only; 
            end
        end
        if SimLys_gteq_expLys == 4
            break
        end
        k = k+1; free_LysError = 5*k;
    end
    
    if (data.protBins_allTimes(end) == data.protBins(end)) && (length(data.protBins) == 2)
        Lys_P_init      = ones(1, 1+size(data.SILAC_data_allTimes,2));
        [t,Lys_P_simu]	= ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final,SILAC_food_impurity),data.t,Lys_P_init);
        Prot_error  	= (nansum((data.SILAC_data_allTimes(:,1:end) - Lys_P_simu(:,2:end)).^2));
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data,  [], SILAC_food_impurity, params.setting);
        Ci              = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
     	Glob_fit_Prot   = [Opts(2:end-1)',((Ci(2:end-1,2)-Ci(2:end-1,1))/2),Prot_error'];
        
    else    
        Glob_fit_Prot = []; 
        for i = 1:length(data.protBins)-1
            fprintf('\nOptimizing proteins-%d to %d, out of %d  with setting-%d\n',data.protBins(i)+1,data.protBins(i+1), data.protBins(end),params.setting	)
            number_param    = 1 + size(data.SILAC_data(:,data.protBins(i)+1:data.protBins(i+1)),2);%
            lb              = zeros(1,number_param); ub = [];
            Ansatz          = [data.Ansatz_(1, 2:number_param) Concen_ratio]; 
          	data.SILAC_data_temp    = [PredLys data.SILAC_data(:,data.protBins(i)+1:data.protBins(i+1))];
        	data.Conc_temp          = [0];
          	if params.opti_algo == 2
                problem     = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, SILAC_food_impurity, params.setting+1, gama_Lys_best),'lb',lb,'ub',ub,'options',opts3);
            else
                problem     = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO_fixLys(param, data.t,  data.SILAC_data_temp,data.Conc_temp, SILAC_food_impurity, params.setting+1, gama_Lys_best),'lb',lb,'ub',ub,'options',opts);
            end
            Ansatz          = [data.Ansatz_(2:end, 2:number_param) repelem(Concen_ratio, size(data.Ansatz_,1)-1)'];
            custpts = CustomStartPointSet(Ansatz);
            [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
            Error_proteins = 100; Gfit_final = []; 
            for j =1 : length(solutions)
                Lys_P_init      = ones(1,number_param);
             	[t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,solutions(j).X, gama_Lys_best,SILAC_food_impurity),data.t,Lys_P_init);
            	Error_proteins_ = sum(nansum((data.SILAC_data_temp(:,2:end) - Lys_P(:,2:end)).^2));
                if Error_proteins_ < Error_proteins
                    Error_proteins  = Error_proteins_; 
                    Gfit_final      = solutions(j).X;
                end
            end
            [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, SILAC_food_impurity, params.setting+1, gama_Lys_best);
            Ci              = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
            Lys_P_init      = ones(1,number_param);
            [t,Lys_P_simu]	= ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts(1:end-1) Opts(end)], gama_Lys_best, SILAC_food_impurity),data.t,Lys_P_init);
            Prot_error      = (nansum((data.SILAC_data_temp(:,2:end) - Lys_P_simu(:,2:end)).^2));
            Glob_fit_Prot(data.protBins(i)+1:data.protBins(i+1),:) = [Opts(1:end-1)',((Ci(1:end-1,2)-Ci(1:end-1,1))/2),Prot_error'];
        end
    end
    Gfit_tab	 = array2table( [data.SILAC_data', Glob_fit_Prot],'VariableNames',{data.timepoints{1:end}, 'HalfLife_in_days', 'Confidence_Interval', 'residual_error' });
 	writetable(struct2table(params),outfile, 'Sheet','parameter_file');
  	writetable([data.ProtInfo,Gfit_tab],outfile, 'Sheet','results');
end
    

if (params.setting	 == 2)
	Glob_fit_Prot = [];     
    for i = 1:length(data.protBins)-1
	    fprintf('\nOptimizing proteins-%d to %d, out of %d  with setting-%d \n',data.protBins(i)+1,data.protBins(i+1), data.protBins(end),params.setting	)
	    number_param	= 2 + size(data.SILAC_data(:,data.protBins(i)+1:data.protBins(i+1)),2);%
	    lb	      = zeros(1,number_param); ub = [];
	    Ansatz	  = data.Ansatz_(1, 1:number_param); 
		data.SILAC_data_temp    = [data.LysRatio' data.SILAC_data(:,data.protBins(i)+1:data.protBins(i+1))];
        if params.opti_algo == 2
            problem             = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,[], SILAC_food_impurity, params.setting),'lb',lb,'ub',ub,'options',opts3);
        else
            problem             = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO_fmincon(param, data.t,  data.SILAC_data_temp,[], SILAC_food_impurity, params.setting),'lb',lb,'ub',ub,'options',opts);
        end
	    Ansatz	= data.Ansatz_(2:end, 1:number_param);
	    custpts	= CustomStartPointSet(Ansatz);
	    [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
        Error_free_Lys = 100;
        for j =1 : length(solutions)
            Lys_P_init      = ones(1,number_param-1);
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X,SILAC_food_impurity),data.t,Lys_P_init);
            Error_free_Lys_ = sum((Lys_P(:,1) - data.LysRatio').^2);
            if (Error_free_Lys_ < Error_free_Lys )
                Gfit_final      = solutions(j).X;
                Error_free_Lys  = Error_free_Lys_;
            end
        end
           
        Lys_P_init              = ones(1,size(data.SILAC_data_temp,2));
        [Opts,~,resid,~,~,~,J]  = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data_temp,  [], SILAC_food_impurity, params.setting	);
        Ci                      = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        [t,Lys_P_simu]          = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,[log(2)./Opts], SILAC_food_impurity),data.t,Lys_P_init);
    	Prot_error              = (nansum((data.SILAC_data_temp(:,2:end) - Lys_P_simu(:,2:end)).^2));
    	Glob_fit_Prot(data.protBins(i)+1:data.protBins(i+1),:) = [Opts(2:end-1)',((Ci(2:end-1,2)-Ci(2:end-1,1))/2),Prot_error'];
    end
    %save the calculated half-lives to output file
 	Gfit_tab            = array2table([ data.SILAC_data',Glob_fit_Prot],'VariableNames',{data.timepoints{1:end}, 'HalfLife_in_days', 'Confidence_Interval', 'residual_error' });
    writetable(struct2table(params),outfile, 'Sheet','parameter_file');
    writetable([data.ProtInfo, Gfit_tab],outfile, 'Sheet','results');
end    

if (params.setting == 3)
	Glob_fit_Prot = []; 
    for i = 1:length(data.protBins)-1
	    fprintf('\nOptimizing proteins-%d to %d, out of %d  with setting-%d \n',data.protBins(i)+1,data.protBins(i+1), data.protBins(end),params.setting)
	    number_param	= 2 + size(data.SILAC_data(:,data.protBins(i)+1:data.protBins(i+1)),2);%
	    lb              = zeros(1,number_param); ub = [];
	    Ansatz          = data.Ansatz_(1, 1:number_param); 

		data.SILAC_data_temp	= [data.LysRatio' data.SILAC_data(:,data.protBins(i)+1:data.protBins(i+1))];
		data.Conc_temp          = [data.freeLysConc, data.ProtConc(data.protBins(i)+1:data.protBins(i+1)), data.TotalLysConc-sum(data.ProtConc(data.protBins(i)+1:data.protBins(i+1)))];
        if params.opti_algo == 2
            problem     = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, SILAC_food_impurity, params.setting),'lb',lb,'ub',ub,'options',opts3);
        else
            problem     = createOptimProblem('fmincon','x0',Ansatz,'objective',@(param)GO_fmincon(param, data.t,  data.SILAC_data_temp,data.Conc_temp, SILAC_food_impurity, params.setting),'lb',lb,'ub',ub,'options',opts);
        end
	    Ansatz          = data.Ansatz_(2:end, 1:number_param);
	    custpts         = CustomStartPointSet(Ansatz);
	    [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
        Error_free_Lys = 100;
        for j =1 : length(solutions)
            Lys_P_init      = ones(1,number_param);
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,solutions(j).X,data.Conc_temp, SILAC_food_impurity),data.t,Lys_P_init);
            Error_free_Lys_ = sum((Lys_P(:,1) - data.LysRatio').^2);
            if (Error_free_Lys_ < Error_free_Lys )
                Gfit_final      = solutions(j).X;
                Error_free_Lys  = Error_free_Lys_;
            end
        end
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, SILAC_food_impurity, params.setting);
        Ci              = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
    	Lys_P_init      = ones(1,1+size(data.SILAC_data_temp,2));
        [t,Lys_P_simu]	= ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,[log(2)./Opts],data.Conc_temp, SILAC_food_impurity),data.t,Lys_P_init);
      	Prot_error      = (nansum((data.SILAC_data_temp(:,2:end) - Lys_P_simu(:,2:end-1)).^2));
    	Glob_fit_Prot(data.protBins(i)+1:data.protBins(i+1),:) = [Opts(2:end-1)',((Ci(2:end-1,2)-Ci(2:end-1,1))/2),Prot_error'];
    end
    %save the calculated half-lives to output file
    synthesis_rate              = (data.ProtConc)'.* (log(2)./Glob_fit_Prot(:,1));
	perce_Frac_synthesis_rate   = (synthesis_rate./data.ProtConc')*100;
	Gfit_tab                    = array2table([ data.SILAC_data', data.ProtConc', Glob_fit_Prot,  synthesis_rate, perce_Frac_synthesis_rate],'VariableNames',{data.timepoints{1:end} 'protein_conc', 'HalfLife_in_days', 'Confidence_Interval', 'residual_error', 'synthesis_rate', 'Perce_Frac_synthesis_rate'});
    writetable(struct2table(params),outfile, 'Sheet','parameter_file');
    writetable([data.ProtInfo, Gfit_tab],outfile, 'Sheet','results');
end

end

