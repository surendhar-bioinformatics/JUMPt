function [OptiResStep1] = gama_Lys(data,parameters,Setting)
%step-1 ;Optimization to get Lys degredation value
%if the proteins are < 100 the half-lives (including CI) are exported to outfile

opts    = optimset('FinDiffRelStep',1e-3, 'MaxIter', 1000);
opts2	= optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',1e-3, 'MaxIterations', 2, 'MaxFunctionEvaluations', 2, 'Display','off');%,
ms      = MultiStart('UseParallel','always','Display','off');
Error_free_Lys =  100; OptiResStep1.gama_Lys_best =  100; min_free_lys_error = 100; SimLys_gteq_expLys =0; free_LysError = 0; k = 0; OptiResStep1.free_Lys_fitting = []; OptiResStep1.gama_Lys_best =  100;  min_Error_proteins_only = 100; min_Error_prot= 100;  
while ((Setting ==3) || (Setting ==2)) && (free_LysError <= 10) %&& (data.se2(end) ~= data.se(end)) ) %|| ((free_LysError < 10) && (data.se2(end) == data.se(end)) )  %&& (OptiResStep1.gama_Lys_best ==  100)
    for i = 1:length(data.se2)-1
        fprintf('\nStep-1; Optimizing proteins-%d to %d, out of %d  with setting-%d \n',data.se2(i)+1,data.se2(i+1), data.se2(end),Setting)
        number_param	= 2 + size(data.SILAC_data2(:,data.se2(i)+1:data.se2(i+1)),2);%
        lb          = zeros(1,number_param); ub = [];
        Ansatz      = data.Ansatz_(1, 1:number_param); 
        fast_data       = data.Lys(2:end)';

        if Setting == 3 
            data.SILAC_data_temp    = [data.Lys'-data.Lys'*free_LysError/100 data.SILAC_data2(:,data.se2(i)+1:data.se2(i+1))];
            data.Conc_temp      = [data.LysConc, data.ProtConc2(data.se2(i)+1:data.se2(i+1)), parameters.tot_Lys-sum(data.ProtConc2(data.se2(i)+1:data.se2(i+1)))];
            problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting),'lb',lb,'ub',ub,'options',opts);
        elseif Setting == 2 
            data.SILAC_data_temp    = [data.Lys'-data.Lys'*free_LysError/100 data.SILAC_data2(:,data.se2(i)+1:data.se2(i+1))];
            data.Conc_temp      = [0];
            problem         = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting),'lb',lb,'ub',ub,'options',opts);
        end

        Ansatz	= data.Ansatz_(2:end, 1:number_param);
        custpts	= CustomStartPointSet(Ansatz);
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
            if ((Setting ==3) && (Error_free_Lys_ < Error_free_Lys )) || ((Setting ==2) && (Error_free_Lys_ < Error_free_Lys ) && (number_param <= 12)) || ((Setting ==2) && (number_param > 12) && ((length(temp(Lys_P(2:end,1) <= fast_data)) > SimLys_gteq_expLys) || ((length(temp(Lys_P(2:end,1) <= fast_data)) == SimLys_gteq_expLys ) && (Error_free_Lys_ < Error_free_Lys ))))
            %if (Error_free_Lys_ < Error_free_Lys )
                SimLys_gteq_expLys = length(temp(Lys_P(2:end,1) <= fast_data));
                OptiResStep1.free_LysError = free_LysError;
                Error_free_Lys = Error_free_Lys_; 
                OptiResStep1.gama_Lys_best =  solutions(j).X(1);
                OptiResStep1.Concen_ratio = solutions(j).X(end);
                Gfit_final = solutions(j).X;
                if Setting == 3
                    [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,solutions(j).X,data.Conc_temp),data.t_long,Lys_P_init);
                    OptiResStep1.free_Lys_fitting  = [[0; t] [Error_free_Lys; Lys_P_simu(:,1)]];
                elseif Setting == 2
                    [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t_long,Lys_P_init);
                    OptiResStep1.free_Lys_fitting  = [[0; t] [Error_free_Lys; Lys_P_simu(:,1)]];
                end
            elseif Error_free_Lys_ < min_free_lys_error
                min_free_lys_error = Error_free_Lys_; 
                gama_Lys_withMinFreeLysError =  solutions(j).X(1); Concen_ratio_withMinFreeLysError = solutions(j).X(end); 
                Gfit_final_withMinFreeLysError = solutions(j).X;
            end
        end
     end
    if SimLys_gteq_expLys == 4
    	break
    end
    k = k+1; free_LysError = 5*k;
end

while (Setting == 1) && (free_LysError <= 20 ) %  
    
    for i = 1:length(data.se2)-1
        fprintf('\nStep-1; Optimizing proteins-%d to %d, out of %d  with setting-%d \n',data.se2(i)+1,data.se2(i+1), data.se2(end),Setting)
        number_param	= 2 + size(data.SILAC_data2(:,data.se2(i)+1:data.se2(i+1)),2);%
        lb              = zeros(1,number_param); ub = [];
        Ansatz          = data.Ansatz_(1, 1:number_param); 

    	data.SILAC_data_temp    = [data.SILAC_data2(:,data.se2(i)+1:data.se2(i+1))];
      	fast_data               = [min(data.SILAC_data2(2,:)) min(data.SILAC_data2(3,:)) min(data.SILAC_data2(4,:)) min(data.SILAC_data2(5,:))]';
    	data.SILAC_data_temp    = [[1; fast_data - fast_data*(free_LysError/100)] data.SILAC_data_temp];
     	data.Conc_temp          = 0;
        problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting+1),'lb',lb,'ub',ub,'options',opts);
        
        Ansatz	= data.Ansatz_(2:end, 1:number_param);
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,~,~,~,solutions]  = run(ms,problem,custpts);
        for j =1 : length(solutions)
            Lys_P_init      = ones(1,number_param-1);
            [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t,Lys_P_init);
         	Error_proteins_only = sum(nansum((data.SILAC_data_temp(:,2:end) - Lys_P(:,2:end)).^2));

            temp = [2; 2; 2; 2];
            if (length(temp(Lys_P(2:end,1) <= fast_data)) > SimLys_gteq_expLys) || ((length(temp(Lys_P(2:end,1) <= fast_data)) == SimLys_gteq_expLys ) && (Error_proteins_only < min_Error_prot )) 
                OptiResStep1.PredLys = Lys_P(:,1);
                SimLys_gteq_expLys = length(temp(Lys_P(2:end,1) <= fast_data));
                min_Error_prot = Error_proteins_only;
                OptiResStep1.free_LysError = free_LysError;
                OptiResStep1.gama_Lys_best =  solutions(j).X(1);
                OptiResStep1.Concen_ratio = solutions(j).X(end);
                Gfit_final = solutions(j).X;
                [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,solutions(j).X),data.t_long,Lys_P_init);
                OptiResStep1.free_Lys_fitting  = [t Lys_P_simu(:,1)];
            elseif Error_proteins_only < min_Error_proteins_only
                min_Error_proteins_only = Error_proteins_only; 
                gama_Lys_withMinFreeLysError =  solutions(j).X(1); Concen_ratio_withMinFreeLysError = solutions(j).X(end); 
                Gfit_final_withMinFreeLysError = solutions(j).X;
            end
        end
    end
    if SimLys_gteq_expLys == 4
        break
    end
    k = k+1; free_LysError = 5*k;
end
if OptiResStep1.gama_Lys_best == 100 
	OptiResStep1.gama_Lys_best = gama_Lys_withMinFreeLysError;
	OptiResStep1.Concen_ratio  = Concen_ratio_withMinFreeLysError;
	Gfit_final    = Gfit_final_withMinFreeLysError;
end

% print to file if there is only one bin to calculate
OptiResStep1.prot_deg_rate = log(2)./Gfit_final(2:end-1)';
if (data.se2(end) == data.se(end)) && (length(data.se2) == 2)
    if Setting == 3
        data.Conc_temp      	= [data.LysConc, data.ProtConc2, parameters.tot_Lys-sum(data.ProtConc2)];
        Lys_P_init              = ones(1,2+size(data.SILAC_data2,2));
        [Opts,~,resid,~,~,~,J]  = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting);
        Ci                      = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        [t,Lys_P_simu]          = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,Gfit_final,data.Conc_temp),data.t,Lys_P_init);
        Prot_error              = (nansum((data.SILAC_data2(:,1:end) - Lys_P_simu(:,2:end-1)).^2));
        Glob_fit_Prot           = [Opts(2:end-1)',(((Ci(2:end-1,2)-Ci(2:end-1,1))/2)/2),Prot_error'];
        synthesis_rate          = (data.ProtConc)'.* (log(2)./Opts(2:end-1))';
        perce_Frac_synthesis_rate = (synthesis_rate./data.ProtConc')*100;
        Gfit_tab                = array2table([ data.SILAC_data', Glob_fit_Prot,  synthesis_rate, perce_Frac_synthesis_rate],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI', 'residual_error', 'synthesis_rate', 'Perce_Frac_synthesis_rate'});
        writetable([data.protInfo, Gfit_tab],parameters.out_file, 'Sheet','Setting3');
    elseif Setting == 2              
        data.Conc_temp              = [];
        Lys_P_init                  = ones(1, 1+size(data.SILAC_data2,2));
        [t,Lys_P_simu]              = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final),data.t,Lys_P_init);
        Prot_error                  = (nansum((data.SILAC_data2(:,1:end) - Lys_P_simu(:,2:end)).^2));
        [Opts,~,resid,~,~,~,J]  = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting);
        Ci                          = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        Glob_fit_Prot               = [Opts(2:end-1)',(((Ci(2:end-1,2)-Ci(2:end-1,1))/2)/2),Prot_error'];
        Gfit_tab                    = array2table([data.SILAC_data',  Glob_fit_Prot],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI', 'residual_error' });
        writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting2');
    elseif Setting == 1
        data.Conc_temp              = [];
        Lys_P_init                  = ones(1, 1+size(data.SILAC_data2,2));
        [t,Lys_P_simu]              = ode15s(@(t,Lys_P)PT_ODE_Ratio(t,Lys_P,Gfit_final),data.t,Lys_P_init);
        Prot_error                  = (nansum((data.SILAC_data2(:,1:end) - Lys_P_simu(:,2:end)).^2));
        [Opts,~,resid,~,~,~,J]  = lsqnonlin(@GO_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data,  data.Conc_temp, Setting);
        Ci                          = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
     	Glob_fit_Prot               = [Opts(2:end-1)',(((Ci(2:end-1,2)-Ci(2:end-1,1))/2)/2),Prot_error'];
        Gfit_tab                    = array2table([data.SILAC_data',  Glob_fit_Prot],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI', 'residual_error' });
        writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting1');
    end
end
end

