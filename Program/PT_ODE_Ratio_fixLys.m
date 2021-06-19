function output1 = PT_ODE_Ratio_fixLys(t,y0, param, gama_Lys, SILAC_food_impurity )
	y1 = y0(1);
	y2 = transpose(y0(2:end)); % row vector

    % define the odes
	dy1 = gama_Lys*(SILAC_food_impurity-y1) + sum(param(1 : end-1).*param(end).*(y2-y1)); %free-Lys
	dy2 = param(1 : end-1) .*(y1 - y2); %proteins
    output1 = [dy1;transpose(dy2)];
end

