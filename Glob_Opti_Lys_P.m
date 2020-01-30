function LSE_Lys_Prot = Glob_Opti_Lys_P(param, tspan, Lys_ratio, P_ratio, Lys_conc, EtaP)
	Lys_P_init = ones(1,1+length(EtaP));
	[t,Lys_P] = ode15s(@(t,y0)PT_ODE(t,y0,param, Lys_conc,EtaP),tspan,Lys_P_init);
    
    LSE_Lys       = sum(nansum((Lys_ratio - Lys_P(:,1)).^2));
    LSE_Lys_Prot  = LSE_Lys + sum(nansum((P_ratio - Lys_P(:,2:end-1)).^2));
    
end