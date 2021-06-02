function LSE_Lys_Prot = GO_fixLys_lsqnonlin2(param, tspan, Lys_P_ratio, EtaP,setting,gama_Lys)
    if setting == 3
        Lys_P_init = ones(1,length(EtaP));
        [t,Lys_P] = ode15s(@(t,y0)PT_ODE_fixLys(t,y0,log(2)./param, EtaP,  gama_Lys),tspan,Lys_P_init );
        LSE_Lys_Prot  = [Lys_P_ratio - Lys_P(:,1:end-1) [0; 0; 0; 0; 0]]; 
        LSE_Lys_Prot = abs(LSE_Lys_Prot); is = isfinite(LSE_Lys_Prot) ;LSE_Lys_Prot(~is) = nanmean(LSE_Lys_Prot(is));
        %LSE_Lys_Prot(isnan(LSE_Lys_Prot))= [];
    end
    if setting == 2
        Lys_P_init = ones(1,size(Lys_P_ratio,2));
        [t,Lys_P] = ode15s(@(t,y0)PT_ODE_Ratio_fixLys(t,y0,[log(2)./param(1:end-1), param(end)],gama_Lys ),tspan,Lys_P_init);
        LSE_Lys_Prot  = [Lys_P_ratio - Lys_P(:,1:end)]; 
        LSE_Lys_Prot = abs(LSE_Lys_Prot); is = isfinite(LSE_Lys_Prot) ;LSE_Lys_Prot(~is) = nanmean(LSE_Lys_Prot(is));
        %LSE_Lys_Prot(isnan(LSE_Lys_Prot))= [];
    end
    
    if setting == 1
        Lys_P_init = ones(1,1+size(Lys_P_ratio,2));
        [t,Lys_P] = ode15s(@(t,y0)PT_ODE_Ratio_fixLys(t,y0,[log(2)./param(1:end-1), param(end)], gama_Lys),tspan,Lys_P_init);
        LSE_Lys_Prot  = [Lys_P_ratio - Lys_P(:,2:end)]; 
        LSE_Lys_Prot = abs(LSE_Lys_Prot); is = isfinite(LSE_Lys_Prot) ;LSE_Lys_Prot(~is) = nanmean(LSE_Lys_Prot(is));
        %LSE_Lys_Prot(isnan(LSE_Lys_Prot))= [];
    end
        
end
