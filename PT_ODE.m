function output1 = PT_ODE(t,y0, param, Lys_conc, EtaP)
	y1 = y0(1);
	y2 = transpose(y0(2:end)); % row vector
    gama_Lys = param(1); 
	gama_P = param(2 : end);
    
	% define the odes
	dy1 = gama_Lys*(0.05-y1) + sum(gama_P.*(EtaP/Lys_conc).*(y2-y1));
	dy2 = gama_P .*(y1 - y2); %proteins
    output1 = [dy1;transpose(dy2)];
end

