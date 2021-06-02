function output1 = PT_ODE_Ratio_fixLys(t,y0, param, gama_Lys)
	y1 = y0(1);
	y2 = transpose(y0(2:end)); % row vector
    %gama_Lys = param(1); 
    gama_P = param(1 : end-1);
    Ratio = param(end);
    %keyboard
	% define the odes
	dy1 = gama_Lys*(0.01-y1) + sum(gama_P.*Ratio.*(y2-y1));
	dy2 = gama_P .*(y1 - y2); %proteins
    output1 = [dy1;transpose(dy2)];
end

