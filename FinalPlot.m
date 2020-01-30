close all 
clearvars

%% read the original data
input_file  = 'test_data_5P.xlsx';
P_data      = readtable(input_file,'Sheet','min3data_orig');%,'Range','1:12',
tspan       = table2array(P_data(1,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'}));
Lys_ratio   = transpose(table2array(P_data(2,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'})));
P_ratio     = transpose(table2array(P_data(3:end,{'data_1' 'data_2' 'data_3' 'data_4' 'data_5'})));
T_array     = transpose(table2array(P_data(3:end,{'T_Day0' 'T_Day4' 'T_Day8' 'T_Day16' 'T_Day32'})));
EtaP        = transpose(table2array(P_data(3:end,'Lys_Conc_microM') ));

Lys_conc    = 206;%microM 
sum_EtaP    = 183190; %micoM
len_P_data  = (size(P_ratio,2)); 
tspan_long  = linspace(0, 32,3201);
tspan_medium = linspace(0, 32,321);
%% Make final plots
gamma_Lys_Prot = readtable('Results\Glob_fit_res.xlsx','Sheet','Glob_fit');
gamma_GlobOpti = transpose((table2array(gamma_Lys_Prot(:,'gamma') )));
EtaP_temp       = [EtaP, sum_EtaP-sum(EtaP)];
%Simulate the data using optimized parameters
Lys_P_init      = ones(1,len_P_data+2);
[t,Lys_P_simu]	= ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,gamma_GlobOpti, Lys_conc, EtaP_temp),tspan,Lys_P_init);

LSE_Lys         = sum(nansum((Lys_ratio(:,1)- Lys_P_simu(:,1)).^2));
LSE_Lys_Prot    = sum(nansum((P_ratio - Lys_P_simu(:,2:end-1)).^2));
[t,Lys_P_simu]  = ode15s(@(t,Lys_P)PT_ODE(t,Lys_P,gamma_GlobOpti, Lys_conc, EtaP_temp),tspan_medium,Lys_P_init);

Prot_name = [];
Prot_name{1}   = 'Lys';
for k = 1:1:len_P_data
    Prot_name{k+1} = sprintf('Protein_%d',k);
end
Prot_name{k+2}     = 'UIP';
f_Lys_Protein_tab = array2table(Lys_P_simu,'VariableNames',Prot_name);
writetable(f_Lys_Protein_tab,'Results\Glob_fit_res.xlsx','Sheet','traj_Lys_P');

% plot the Fitting results
figure('visible','off'); hold on; box on; grid on; set(gca,'LineWidth',2.5); set(gca,'FontSize',14)
plot(t, Lys_P_simu(:,1),'-r','MarkerSize',8,'MarkerEdgeColor','g','LineWidth',2);
plot(t, Lys_P_simu(:,2:end),'-b','MarkerSize',8,'LineWidth',1.5);
for i=1:len_P_data
    plot((T_array((~isnan(T_array(1:5,i))),i)), (P_ratio((~isnan(P_ratio(1:5,i))),i)),'sb','MarkerSize',8,'LineWidth',1.5);
end
plot(t,Lys_P_simu(:,1),'-r','MarkerSize',8,'MarkerEdgeColor','r','LineWidth',2);
plot(tspan, Lys_ratio,'rs','MarkerSize',8,'MarkerEdgeColor','r','LineWidth',2);
title('Global Fit of Lys and all proteins','FontSize',14);
xlabel('Time (days)', 'FontSize',14); ylabel('Fraction of  ''light''  Lys ','FontSize',14);
hlegend= legend('A_L- Fitting','P_L- Fitting'); set(hlegend, 'FontSize',14, 'Box','on'); ylim([0 1.03]);hold off
print('Results/GlobalFit','-dpng','-r400'); hold off
