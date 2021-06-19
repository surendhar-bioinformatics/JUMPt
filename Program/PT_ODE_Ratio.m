function output1 = PT_ODE_Ratio(t,y0, param, SILAC_food_impurity)
% ODE for free-Lys and proteins
% Average ratio between free-Lys and protein bound Lys

	y1 = y0(1);
	y2 = transpose(y0(2:end)); % row vector
    
	% define the odes
	dy1 = param(1)*(SILAC_food_impurity-y1) + sum(param(2 : end-1).*param(end).*(y2-y1)); %free-Lys
	dy2 = param(2 : end-1) .*(y1 - y2); %proteins
    output1 = [dy1;transpose(dy2)];
end

