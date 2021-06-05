function [] = Prot_turnover_v1(data,parameters)
%program to study the protein turnover kinetics. version 1.0.0
 
%read SILAC data 
size_Protein_data = size(data.SILAC_data,2);
fun_deg = @(b,tspan) exp(-b*tspan); 
for i =1:size_Protein_data
    Pdata_temp = data.SILAC_data(:,i);
	non_NaN_index_of_Pdata_temp = find(~isnan(Pdata_temp));
	Pdata_temp = Pdata_temp(non_NaN_index_of_Pdata_temp);
    rng default % For reproducibility
	ExpDegRate_temp = fminsearch(@(b) norm(Pdata_temp - fun_deg(b,data.t(non_NaN_index_of_Pdata_temp))), rand(1,1)) ;
	apparent_half_lives(i) = ExpDegRate_temp;
end
apparent_half_lives = transpose(log(2)./apparent_half_lives);
[~, idx] = sort(apparent_half_lives);
data.SILAC_data = data.SILAC_data(:,idx); 
data.ProtConc = data.ProtConc(:,idx);
data.protInfo = data.protInfo(idx, :); 

data.SILAC_data2 = data.SILAC_data(:, (sum(~isnan(data.SILAC_data)) > 4));
data.protInfo2 = data.protInfo((sum(~isnan(data.SILAC_data)) > 4),:);
data.ProtConc2 = data.ProtConc(:, (sum(~isnan(data.SILAC_data)) > 4));

% define Bin_size and Bin#
if parameters.Num_prot_to_fit < parameters.prot_per_fit
    Bin_size = parameters.Num_prot_to_fit;
else
    Bin_size = parameters.prot_per_fit;
end

Bins = ceil(size_Protein_data/Bin_size);
indices = ([0:Bin_size-1]*Bins)+1;
indices_0 = indices;
for i = 1:Bins-1
    indices = [indices (indices_0+i)];
end
indices(indices > size_Protein_data) = [];
data.SILAC_data = data.SILAC_data(:,indices); 
data.ProtConc = data.ProtConc(:,indices);
data.protInfo = data.protInfo(indices, :); 

data.SILAC_data = data.SILAC_data(:, 1:parameters.Num_prot_to_fit); 
data.ProtConc = data.ProtConc(:, 1:parameters.Num_prot_to_fit);
data.protInfo = data.protInfo(1:parameters.Num_prot_to_fit,:); 
size_Protein_data = size(data.SILAC_data,2);
%creating an array with bin size
se = [0];
for i = 1:floor(size_Protein_data/parameters.prot_per_fit)
    se(i+1) = i*parameters.prot_per_fit;
end
if (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) > 50)
    se(end+1) = size_Protein_data; 
elseif (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) < 51) && (length(se) >1)
    se(end) = size_Protein_data; 
elseif (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) < 51) && (length(se) == 1)
    se(end+1) = size_Protein_data; 
end

%% create araay and bin for selcted protein with all time points 
size_Protein_data2 = size(data.SILAC_data2,2);
% define Bin_size and Bin#
if parameters.Num_prot_to_fit < parameters.prot_per_fit
    Bin_size = parameters.Num_prot_to_fit;
else
    Bin_size = parameters.prot_per_fit;
end

Bins = ceil(size_Protein_data2/Bin_size);
indices = ([0:Bin_size-1]*Bins)+1;
indices_0 = indices;
for i = 1:Bins-1
    indices = [indices (indices_0+i)];
end
indices(indices > size_Protein_data2) = [];
data.SILAC_data2 = data.SILAC_data2(:,indices); 
data.ProtConc2 = data.ProtConc2(:,indices);
data.protInfo2 = data.protInfo2(indices, :);

if parameters.Num_prot_to_fit < size_Protein_data2
    data.SILAC_data2 = data.SILAC_data2(:, 1:parameters.Num_prot_to_fit); 
    data.ProtConc2 = data.ProtConc2(:, 1:parameters.Num_prot_to_fit);
    data.protInfo2 = data.protInfo2(1:parameters.Num_prot_to_fit, :); 
end
size_Protein_data2 = size(data.SILAC_data2,2);
%creating an array with bin size
se2 = [0];
for i = 1:floor(size_Protein_data2/parameters.prot_per_fit)
    se2(i+1) = i*parameters.prot_per_fit;
end
if (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) > 50)
    se2(end+1) = size_Protein_data2; 
elseif (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) < 51) && (length(se2) >1)
    se2(end) = size_Protein_data2; 
elseif (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) < 51) && (length(se2) == 1)
    se2(end+1) = size_Protein_data2; 
end

if length(se2) > 5
    se2 = se2(1:5);
end 

rng default % For reproducibility
Ansatz_          = rand(11,size_Protein_data+2);


%% Calculate the half-lives with different settings
if parameters.setting3 == 1
    [gama_Lys_best,Concen_ratio, free_Lys_fitting] = gama_Lys(se2,se, data,parameters,Ansatz_,3);
    if se(end) >  se2(end)
        Glob_fit_Prot_setting3 = gama_Prot(se,data,parameters,Ansatz_, gama_Lys_best, Concen_ratio, free_Lys_fitting, 3);
    end
end
if parameters.setting2 == 1
    [gama_Lys_best,Concen_ratio, free_Lys_fitting] = gama_Lys(se2,se, data,parameters,Ansatz_,2);
    if se(end) >  se2(end)
        Glob_fit_Prot_setting2 = gama_Prot(se,data,parameters,Ansatz_, gama_Lys_best, Concen_ratio, free_Lys_fitting, 2);
    end
end
if parameters.setting1 == 1
    [gama_Lys_best,Concen_ratio, free_Lys_fitting] = gama_Lys(se2, se, data,parameters,Ansatz_,1);
    if se(end) >  se2(end)
        Glob_fit_Prot_setting1 = gama_Prot(se,data,parameters,Ansatz_, gama_Lys_best, Concen_ratio, free_Lys_fitting, 1);
    end
end
fprintf('\n *****  Completed exporting half-live to excel file ; Now the Program is complete *******\n\n')
end