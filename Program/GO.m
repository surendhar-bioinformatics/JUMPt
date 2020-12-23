function LSE_Lys_P = GO(param, tspan, Lys_P_ratio,EtaP)
    Lys_P_init = ones(1,length(EtaP));
	[t,Lys_P] = ode15s(@(t,y0)PT_ODE(t,y0,param,EtaP),tspan,Lys_P_init);
    LSE_Lys_P  = (sum(nansum((Lys_P_ratio - Lys_P(:,1:end-1)).^2)));
end