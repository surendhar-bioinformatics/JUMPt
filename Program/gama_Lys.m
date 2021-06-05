function [gama_Lys_best,Concen_ratio, free_Lys_fitting] = gama_Lys(se2,se, data,parameters,Ansatz_,Setting)
%step-1 ;Optimization to get Lys degredation value
%if the proteins are < 100 the half-lives (including CI) are exported to outfile

opts = optimset('FinDiffRelStep',1e-3);
opts2           = optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',1e-3,'FiniteDifferenceType','central', 'MaxIterations', 2, 'MaxFunctionEvaluations', 2);%,
ms              = MultiStart('UseParallel','always','Display','iter');
if (Setting ==3) || (Setting ==2)
    Error_free_Lys =  100; gama_Lys_best =  100;  min_free_lys_error = 100; free_LysError = 0; k = 0; free_Lys_fitting = [];  
    while (gama_Lys_best == 100) && (free_LysError <= 10) %  
        for i = 1:length(se2)-1
            fprintf('Step-1; Optimizing proteins-%d to %d, out of %d  with setting-%d \n',se2(i)+1,se2(i+1), se2(end),Setting)
            number_param    = 2 + size(data.SILAC_data2(:,se2(i)+1:se2(i+1)),2);%
            lb              = zeros(1,number_param); 
            ub              = [];
            Ansatz          = Ansatz_(1, 1:number_param); 
            fast_data               = data.Lys(2:end)'+ data.Lys(2:end)'*0.05;

            if Setting == 3 
                data.SILAC_data_temp    = [data.Lys'-data.Lys'*free_LysError/100 data.SILAC_data2(:,se2(i)+1:se2(i+1))];
                data.Conc_temp          = [data.LysConc, data.ProtConc2(se2(i)+1:se2(i+1)), parameters.tot_Lys-sum(data.ProtConc2(se2(i)+1:se2(i+1)))];
                problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting),'lb',lb,'ub',ub,'options',opts);
            elseif Setting == 2 
                data.SILAC_data_temp    = [data.Lys'-data.Lys'*free_LysError/100 data.SILAC_data2(:,se2(i)+1:se2(i+1))];
                data.Conc_temp          = [0];
                problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting),'lb',lb,'ub',ub,'options',opts);
            else
                fprintf('Error: Setting is not defined; please define the setting to calculate the half-lives\n')
            end

            Ansatz          = Ansatz_(2:end, 1:number_param);
            custpts = CustomStartPointSet(Ansatz);
            [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);

            for j =1 : length(solutions)
                if Setting == 3
                    Lys_P_init      = ones(1,number_param);
                    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,solutions(j).X,data.Conc_temp),data.t,Lys_P_init);
                    Error_free_Lys_ = sum((Lys_P(:,1) - data.Lys').^2);
                elseif Setting == 2
                    Lys_P_init      = ones(1,number_param-1);
                    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t,Lys_P_init);
                    Error_free_Lys_ = sum((Lys_P(:,1) - data.Lys').^2);
                end
                temp = [2; 2; 2; 2];
                if (Setting == 3) || (Setting == 2)
                    if (length(temp(Lys_P(2:end,1) <= fast_data)) > 2) && (Error_free_Lys_ < Error_free_Lys ) 
                        Error_free_Lys = Error_free_Lys_; 
                        gama_Lys_best =  solutions(j).X(1);
                        Concen_ratio = solutions(j).X(end);
                        Gfit_final = solutions(j).X;
                        if Setting == 3
                            [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,solutions(j).X,data.Conc_temp),data.t_long,Lys_P_init);
                            free_Lys_fitting  = [[0; t] [Error_free_Lys; Lys_P_simu(:,1)]];
                        elseif Setting == 2
                            [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t_long,Lys_P_init);
                            free_Lys_fitting  = [[0; t] [Error_free_Lys; Lys_P_simu(:,1)]];
                        end
                    elseif Error_free_Lys_ < min_free_lys_error
                        min_free_lys_error = Error_free_Lys_; 
                        gama_Lys_withMinFreeLysError =  solutions(j).X(1); Concen_ratio_withMinFreeLysError = solutions(j).X(end); 
                        Gfit_final_withMinFreeLysError = solutions(j).X;
                    end
                end
            end
        end
        k = k+1; free_LysError = 5*k;
    end
end

if (Setting == 1)
    Error_free_Lys =  100; gama_Lys_best =  100;  min_prot_Lys_error = 0; min_prot_error =100; free_LysError = 0; k = 0; free_Lys_fitting = []; 
    while (free_LysError <= 10) %  
        for i = 1:length(se2)-1
            fprintf('Step-1; Optimizing proteins-%d to %d, out of %d  with setting-%d \n',se2(i)+1,se2(i+1), se2(end),Setting)
            number_param    = 2 + size(data.SILAC_data2(:,se2(i)+1:se2(i+1)),2);%
            lb              = zeros(1,number_param); 
            ub              = [];
            Ansatz          = Ansatz_(1, 1:number_param); 

        	data.SILAC_data_temp    = [data.SILAC_data2(:,se2(i)+1:se2(i+1))];
          	fast_data               = [min(data.SILAC_data2(2,:)) min(data.SILAC_data2(3,:)) min(data.SILAC_data2(4,:)) min(data.SILAC_data2(5,:))]';
          	fast_data               = fast_data - fast_data*(free_LysError/100);
        	data.SILAC_data_temp    = [[1; fast_data] data.SILAC_data_temp];
         	data.Conc_temp          = [0];
            problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting+1),'lb',lb,'ub',ub,'options',opts);
            
            Ansatz          = Ansatz_(2:end, 1:number_param);
            custpts = CustomStartPointSet(Ansatz);
            [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);

            for j =1 : length(solutions)
                Lys_P_init      = ones(1,number_param-1);
                [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t,Lys_P_init);
            	Error_proteins_only = sum(nansum((data.SILAC_data_temp(:,2:end) - Lys_P(:,2:end)).^2));
                if Error_proteins_only < min_prot_error 
                    min_prot_error = Error_proteins_only;
                    Error_proteins_ =  solutions(j).Fval;  
                    gama_Lys_best =  solutions(j).X(1);
                    Concen_ratio = solutions(j).X(end);
                    Gfit_final = solutions(j).X;
                    [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t_long,Lys_P_init);
                    free_Lys_fitting  = [t Lys_P_simu(:,1)];
                elseif solutions(j).Fval < min_prot_Lys_error
                    min_prot_Lys_error = solutions(j).Fval;
                    gama_Lys_withMinFreeLysError =  solutions(j).X(1); Concen_ratio_withMinFreeLysError = solutions(j).X(end); 
                    Gfit_final_withMinFreeLysError = solutions(j).X;
                end
            end
        end
        k = k+1; free_LysError = 5*k;
    end
end
if gama_Lys_best == 100 
	gama_Lys_best = gama_Lys_withMinFreeLysError;
	Concen_ratio  = Concen_ratio_withMinFreeLysError;
	Gfit_final    = Gfit_final_withMinFreeLysError;
end
if se2(end) == se(end)
    if Setting == 3
        data.Conc_temp	= [data.LysConc, data.ProtConc2, parameters.tot_Lys-sum(data.ProtConc2)];
        Lys_P_init      = ones(1,2+size(data.SILAC_data2,2));
        [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,Gfit_final,data.Conc_temp),data.t_long,Lys_P_init);
        free_Lys_fitting  = [[0; t] [Error_free_Lys; Lys_P_simu(:,1)]];
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,Gfit_final,data.Conc_temp),data.t,Lys_P_init);
        Prot_error       = (nansum((data.SILAC_data2(:,1:end) - Lys_P_simu(:,2:end-1)).^2));
        Glob_fit_Prot = [Opts(2:end-1)',(((Ci(2:end-1,2)-Ci(2:end-1,1))/2)/2),Prot_error'];
        synthesis_rate = (data.ProtConc)'.* (log(2)./Opts(2:end-1))';
        perce_Frac_synthesis_rate = (synthesis_rate./data.ProtConc')*100;
        Gfit_tab    = array2table([ data.SILAC_data', Glob_fit_Prot,  synthesis_rate, perce_Frac_synthesis_rate],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI', 'residual_error', 'synthesis_rate', 'Perce_Frac_synthesis_rate'});
        writetable([data.protInfo, Gfit_tab],parameters.out_file, 'Sheet','Setting3');
    elseif Setting == 2              
        data.Conc_temp	= [];
        Lys_P_init      = ones(1, 1+size(data.SILAC_data2,2));
        [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final),data.t_long,Lys_P_init);
        free_Lys_fitting  = [[0; t] [Error_free_Lys; Lys_P_simu(:,1)]];
        [t,Lys_P_simu]    = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final),data.t,Lys_P_init);
        Prot_error        = (nansum((data.SILAC_data2(:,1:end) - Lys_P_simu(:,2:end)).^2));
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting);
        Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        Glob_fit_Prot = [Opts(2:end-1)',(((Ci(2:end-1,2)-Ci(2:end-1,1))/2)/2),Prot_error'];
        Gfit_tab    = array2table([data.SILAC_data',  Glob_fit_Prot],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI', 'residual_error' });
        writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting2');
    elseif Setting == 1
        data.Conc_temp	= [];
        Lys_P_init      = ones(1, 1+size(data.SILAC_data2,2));
        [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final),data.t_long,Lys_P_init);
        free_Lys_fitting  = [t Lys_P_simu(:,1)];
        [Opts,~,resid,~,~,~,J] = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data,  data.Conc_temp, Setting);
         Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        [t,Lys_P_simu]    = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final),data.t,Lys_P_init);
        Prot_error        = (nansum((data.SILAC_data2(:,1:end) - Lys_P_simu(:,2:end)).^2));
        Glob_fit_Prot= [Opts(2:end-1)',(((Ci(2:end-1,2)-Ci(2:end-1,1))/2)/2),Prot_error'];
        Gfit_tab    = array2table([data.SILAC_data',  Glob_fit_Prot],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI', 'residual_error' });
        writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting1');
    end
end
end

