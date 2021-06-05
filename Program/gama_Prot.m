function [Glob_fit_Prot] = gama_Prot(se,data,parameters,Ansatz_,gama_Lys_best, Concen_ratio, free_Lys_fitting,Setting)
    %step-1 ;Optimization to get gama_Lys value
    %opts = optimset('FinDiffRelStep',1e-3);
    opts            = optimoptions(@lsqnonlin, 'FiniteDifferenceStepSize',1e-3,'FiniteDifferenceType','central', 'MaxIterations', 50000);%,, 'MaxFunctionEvaluations', 5000
    opts2           = optimoptions(@lsqnonlin,'FiniteDifferenceStepSize',1e-3,'FiniteDifferenceType','central', 'MaxIterations', 2, 'MaxFunctionEvaluations', 2);%,
    ms              = MultiStart('UseParallel','always','Display','iter');
    Glob_fit_Prot = []; Glob_fit_Lys = [];  Glob_fit_UIP = []; res_error_ = []; min_free_lys_error = 100;
    for i = 1:length(se)-1
        fprintf('Step-2; Optimizing proteins-%d to %d, out of %d  with setting-%d\n',se(i)+1,se(i+1), se(end),Setting)
        number_param    = 1 + size(data.SILAC_data(:,se(i)+1:se(i+1)),2);%
        lb              = zeros(1,number_param); 
        ub              = [];
        Ansatz          = [Ansatz_(1, 2:number_param) Concen_ratio]; 
        fast_data               = data.Lys(2:end)'+ data.Lys(2:end)'*0.05;
        if Setting == 3 
            data.SILAC_data_temp    = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
            data.Conc_temp          = [data.LysConc, data.ProtConc(se(i)+1:se(i+1)), parameters.tot_Lys-sum(data.ProtConc(se(i)+1:se(i+1)))];
            problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting, gama_Lys_best),'lb',lb,'ub',ub,'options',opts);
        elseif Setting == 2 
            data.SILAC_data_temp    = [data.Lys' data.SILAC_data(:,se(i)+1:se(i+1))];
            data.Conc_temp          = [0];
            problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting, gama_Lys_best),'lb',lb,'ub',ub,'options',opts);
        elseif Setting == 1
            data.SILAC_data_temp    = [data.SILAC_data(:,se(i)+1:se(i+1))];
            data.Conc_temp          = [0];
            problem                 = createOptimProblem('lsqnonlin','x0',Ansatz,'objective',@(param)GO_fixLys_lsqnonlin(param, data.t,  data.SILAC_data_temp,data.Conc_temp, Setting, gama_Lys_best),'lb',lb,'ub',ub,'options',opts);
        else
            fprintf('Error: Setting is not defined; please define the setting to calculate the half-lives\n')
        end
        Ansatz          = [Ansatz_(2:end, 2:number_param) repelem(Concen_ratio, size(Ansatz_,1)-1)'];
        custpts = CustomStartPointSet(Ansatz);
        [Gfit,fval,exitflag,output,solutions]  = run(ms,problem,custpts);
        if (Setting == 3) || (Setting == 2)
            Error_free_Lys = 100;
            Gfit_final = [];
            for j =1 : length(solutions)
                if Setting == 3
                    Lys_P_init      = ones(1,number_param+1);
                    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_fixLys(t,Lys_P,solutions(j).X,data.Conc_temp,gama_Lys_best),data.t,Lys_P_init);
                    Error_free_Lys_ = sum((Lys_P(:,1) - data.Lys').^2);
                elseif Setting == 2
                    Lys_P_init      = ones(1,number_param);
                    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,solutions(j).X, gama_Lys_best),data.t,Lys_P_init);
                    Error_free_Lys_ = sum((Lys_P(:,1) - data.Lys').^2);
                elseif Setting == 1
                    Lys_P_init      = ones(1,number_param);
                    [t,Lys_P]       = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,solutions(j).X, gama_Lys_best),data.t,Lys_P_init);
                    fast_data        = [min(data.SILAC_data(2,:)) min(data.SILAC_data(3,:)) min(data.SILAC_data(4,:)) min(data.SILAC_data(5,:))]';
                    Error_free_Lys_ = sum((Lys_P(2:end,1) - fast_data).^2);
                end
                temp = [2; 2; 2; 2];
                if (length(temp(Lys_P(2:end,1) <= fast_data)) > 2) && (Error_free_Lys_ < Error_free_Lys ) 
                    Error_free_Lys = Error_free_Lys_; 
                    Gfit_final = solutions(j).X;
                elseif Error_free_Lys_ < min_free_lys_error
                    min_free_lys_error = Error_free_Lys_;
                    Gfit_withMinFreeLysError =  solutions(j).X;
                end
            end
            if length(Gfit_final) < 1 
                Gfit_final = Gfit_withMinFreeLysError;
            end
        else
            Gfit_final = Gfit;
        end
        if Gfit_final(end) == 0; Gfit_final(end)  =1e-32; end;
        
        %calculating the residual error
        if Setting == 3
            [Opts,fval2,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2,[log(2)./Gfit_final],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting, gama_Lys_best);
            Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        
            [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_fixLys(t,Lys_P,[log(2)./Opts],data.Conc_temp, gama_Lys_best),data.t,Lys_P_init);
            Prot_error       = (nansum((data.SILAC_data_temp(:,2:end) - Lys_P_simu(:,2:end-1)).^2));
            free_lys_error = nansum((data.Lys' - Lys_P_simu(:,1)).^2);
            res_error_(i,:) = [sum(Prot_error) free_lys_error];
            Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(1:end-1)',Ci(1:end-1,:),Prot_error'];
            Glob_fit_Lys(i,:) = log(2)/gama_Lys_best;
            Glob_fit_UIP(i,:) = Opts(end);
            [t,Lys_P_simu]       = ode15s(@(t,Lys_P)PT_ODE_fixLys(t,Lys_P,[log(2)./Opts],data.Conc_temp, gama_Lys_best),data.t_long,Lys_P_init);
            free_Lys_fitting = [free_Lys_fitting [free_lys_error; Lys_P_simu(:,1)]];
            prot_fitting(:,se(i)+1:se(i+1)) = [Prot_error; Lys_P_simu(:,2:end-1)];
        elseif Setting == 2
            [Opts,fval2,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting, gama_Lys_best);
            Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
        
            [t,Lys_P_simu]  	= ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts], gama_Lys_best),data.t,Lys_P_init);
            Prot_error          = (nansum((data.SILAC_data_temp(:,2:end) - Lys_P_simu(:,2:end)).^2));
            free_lys_error      = nansum((data.Lys' - Lys_P_simu(:,1)).^2);
            res_error_(i,:)     = [sum(Prot_error) free_lys_error];
            Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(1:end-1)',Ci(1:end-1,:),Prot_error'];
            Glob_fit_Lys(i,:)   = log(2)/gama_Lys_best;
            Glob_fit_ratioConst(i,:) = Opts(end);
            [t,Lys_P_simu]   	= ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts(1:end-1) Opts(end)],gama_Lys_best),data.t_long,Lys_P_init);
            free_Lys_fitting    = [free_Lys_fitting [free_lys_error; Lys_P_simu(:,1)]];
            prot_fitting(:,se(i)+1:se(i+1)) = [Prot_error; Lys_P_simu(:,2:end)];
        elseif Setting == 1
            [Opts,fval2,resid,~,~,~,J] = lsqnonlin(@GO_fixLys_lsqnonlin2,[log(2)./Gfit_final(1:end-1) Gfit_final(end)],lb,ub,opts2, data.t, data.SILAC_data_temp,  data.Conc_temp, Setting, log(2)/gama_Lys_best);
            Ci = nlparci(Opts,resid,'jacobian',J,'alpha',.05);
            Lys_P_init      = ones(1,number_param);
            [t,Lys_P_simu]      = ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts(1:end-1) Opts(end)], gama_Lys_best),data.t,Lys_P_init);
            Prot_error          = (nansum((data.SILAC_data_temp - Lys_P_simu(:,2:end)).^2));
            res_error_(i,:)     = sum(Prot_error);
            Glob_fit_Prot(se(i)+1:se(i+1),:) = [Opts(1:end-1)',Ci(1:end-1,:),Prot_error'];
            Glob_fit_Lys(i,:)   = log(2)/gama_Lys_best;
            Glob_fit_ratioConst(i,:) = Opts(end);
            [t,Lys_P_simu]   	= ode15s(@(t,Lys_P)PT_ODE_Ratio_fixLys(t,Lys_P,[log(2)./Opts(1:end-1) Opts(end)],gama_Lys_best),data.t_long,Lys_P_init);
            free_Lys_fitting    = [free_Lys_fitting, Lys_P_simu(:,1)];
            prot_fitting(:,se(i)+1:se(i+1)) = [Prot_error; Lys_P_simu(:,2:end)];
        end
    end
    if Setting == 3
        synthesis_rate = (data.ProtConc)'.* (log(2)./Glob_fit_Prot(:,1));
        perce_Frac_synthesis_rate = (synthesis_rate./data.ProtConc')*100;
        Gfit_tab    = array2table([ data.SILAC_data', Glob_fit_Prot,  synthesis_rate, perce_Frac_synthesis_rate, data.ProtConc' ],'VariableNames',{'day0','day4','day8','day16','day32', 'HalfLife_in_days', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error', 'synthesis_rate', 'Perce_Frac_synthesis_rate', 'protein_conc'});
        writetable([data.protInfo, Gfit_tab],parameters.out_file, 'Sheet','Setting3');
        Lys_UIP_degRate_tab    = array2table([Glob_fit_Lys Glob_fit_UIP res_error_],'VariableNames',{ 'Lys_degRate','UIP_degRate', 'res_error_prot', 'res_error_Lys' });
        writetable(Lys_UIP_degRate_tab,parameters.out_file, 'Sheet','Setting3_LysDeg');
        writetable(array2table(prot_fitting), parameters.out_file,'Sheet','Setting3_prot_fit');
        writetable(array2table(free_Lys_fitting),parameters.out_file, 'Sheet','Setting3_Lys_fit');
    elseif Setting == 2
        Gfit_tab    = array2table(Glob_fit_Prot,'VariableNames',{'HalfLife_in_days', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error' });
        writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting2');
        Lys_UIP_degRate_tab    = array2table([Glob_fit_Lys Glob_fit_ratioConst res_error_],'VariableNames',{ 'Lys_degRate','ratioConst', 'res_error_prot', 'res_error_Lys' });
        writetable(Lys_UIP_degRate_tab,parameters.out_file, 'Sheet','Setting2_LysDeg');
        writetable(array2table(prot_fitting), parameters.out_file,'Sheet','Setting2_prot_fit');
        writetable(array2table(free_Lys_fitting),parameters.out_file, 'Sheet','Setting2_Lys_fit');
   elseif Setting == 1
        Gfit_tab    = array2table(Glob_fit_Prot,'VariableNames',{'HalfLife_in_days', 'HalfLife_CI_1', 'HalfLife_CI_2', 'residual_error' });
        writetable([data.protInfo,Gfit_tab],parameters.out_file, 'Sheet','Setting1');
        Lys_UIP_degRate_tab    = array2table([Glob_fit_Lys Glob_fit_ratioConst res_error_],'VariableNames',{ 'Lys_degRate','ratioConst', 'res_error_prot' });
        writetable(Lys_UIP_degRate_tab,parameters.out_file, 'Sheet','Setting1_LysDeg');
        writetable(array2table(prot_fitting), parameters.out_file,'Sheet','Setting1_prot_fit');
        writetable(array2table(free_Lys_fitting),parameters.out_file, 'Sheet','Setting1_Lys_fit');
    end
end

