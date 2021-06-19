function LSE_Lys_P = GO_fmincon(param, tspan, Lys_P_ratio,EtaP,SILAC_food_impurity, setting)
    if setting == 3
        Lys_P_init = ones(1,length(EtaP));
        [t,Lys_P] = ode15s(@(t,y0)PT_ODE(t,y0,param,EtaP,SILAC_food_impurity),tspan,Lys_P_init);
        LSE_Lys_P  = (sum(nansum((Lys_P_ratio - Lys_P(:,1:end-1)).^2)));
    end
    if setting == 2
        Lys_P_init = ones(1,size(Lys_P_ratio,2));
        [t,Lys_P] = ode15s(@(t,y0)PT_ODE_Ratio(t,y0,param,SILAC_food_impurity),tspan,Lys_P_init);
        LSE_Lys_P  = (sum(nansum((Lys_P_ratio - Lys_P(:,1:end)).^2)));
    end
    if setting == 1
        Lys_P_init = ones(1,1+size(Lys_P_ratio,2));
        [t,Lys_P] = ode15s(@(t,y0)PT_ODE_Ratio(t,y0,param,SILAC_food_impurity),tspan,Lys_P_init);
        LSE_Lys_P  = (sum(nansum((Lys_P_ratio - Lys_P(:,2:end)).^2)));
    end
end