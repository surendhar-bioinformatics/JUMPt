function output1 = PT_ODE_fixLys(t,y0, gama, EtaP, gama_Lys)
	y1 = y0(1);
	y2 = transpose(y0(2:end)); 
	% define the odes
	dy1 = gama_Lys*(0.01-y1) + sum(gama(1:end).*(EtaP(2:end)/EtaP(1)).*(y2-y1));%Lys
	dy2 = gama(1:end).*(y1 - y2); %proteins
    output1 = [dy1;transpose(dy2)];
end
