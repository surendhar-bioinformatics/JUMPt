function LSE_Lys_Prot = GO_lsqnonlin(param, tspan, Lys_P_ratio, EtaP)
    Lys_P_init = ones(1,length(EtaP));
    %[t,Lys_P] = ode15s(@(t,y0)PT_ODE_halflife(t,y0,param, EtaP),tspan,Lys_P_init );
    [t,Lys_P] = ode15s(@(t,y0)PT_ODE(t,y0,log(2)./param, EtaP),tspan,Lys_P_init );

    LSE_Lys_Prot  = [Lys_P_ratio - Lys_P(:,1:end-1) [0; 0; 0; 0; 0]]; 
    LSE_Lys_Prot(isnan(LSE_Lys_Prot))= [];
end
