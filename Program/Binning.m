function data = Binning(parameters)
%% read the data
All_data    = readtable(parameters.input_file);
mask        = startsWith(All_data.Properties.VariableNames, 'data');
mask3       = startsWith(All_data.Properties.VariableNames, 'Concentration');
mask2       = startsWith(All_data.Name, 'Lys');
data.t      = table2array(All_data(1,mask));  
data.t_long = linspace(0, 32,321);

if any(mask2)
    data.Lys        = table2array(All_data(mask2,mask));% Lys data
    data.protInfo   = (All_data(3:end,1:find(mask,1)-1)); % protein IDs
    data.ProtConc   = table2array(All_data(3:end,mask3))';% Protein concentration
    SILAC_data      = table2cell(All_data(3:end,mask));% protein data
    SILAC_data_NaN  = cellfun(@(x) ~isa(x,'double'),SILAC_data); % check for possible 0x0 char reads
    data.LysConc  = table2array(All_data(2,mask3));
    SILAC_data(SILAC_data_NaN)   = {NaN};% set those to NaN
    data.SILAC_data              = cell2mat(SILAC_data)'; % heavy/light ratios,
else
    data.protInfo  = (All_data(2:end,1:find(mask,1)-1));
    data.ProtConc  = table2array(All_data(2:end,mask3))'; 
    SILAC_data     = table2cell(All_data(2:end,mask));
    SILAC_data_NaN = cellfun(@(x) ~isa(x,'double'),SILAC_data); 
    SILAC_data(SILAC_data_NaN)   = {NaN};
    data.SILAC_data              = cell2mat(SILAC_data)'; 
end

parameters.Num_prot_to_fit  = size(data.SILAC_data,2);

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
data.se = [0];
for i = 1:floor(size_Protein_data/parameters.prot_per_fit)
    data.se(i+1) = i*parameters.prot_per_fit;
end
if (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) >= parameters.prot_per_fit/4)
    data.se(end+1) = size_Protein_data; 
elseif (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) < parameters.prot_per_fit/4) && (length(data.se) >1)
    data.se(end) = size_Protein_data; 
elseif (mod(size_Protein_data,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data,parameters.prot_per_fit) < parameters.prot_per_fit/4) && (length(data.se) == 1)
    data.se(end+1) = size_Protein_data; 
end


%% create araay and bin for selcted protein with all time points 
size_Protein_data2 = size(data.SILAC_data2,2);
%parameters.prot_per_fit = 100;
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
data.se2 = [0];
for i = 1:floor(size_Protein_data2/parameters.prot_per_fit)
    data.se2(i+1) = i*parameters.prot_per_fit;
end
if (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) >= parameters.prot_per_fit/4)
    data.se2(end+1) = size_Protein_data2; 
elseif (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) < parameters.prot_per_fit/4) && (length(data.se2) >1)
    data.se2(end) = size_Protein_data2; 
elseif (mod(size_Protein_data2,parameters.prot_per_fit) > 0)  && (mod(size_Protein_data2,parameters.prot_per_fit) < parameters.prot_per_fit/4) && (length(data.se2) == 1)
    data.se2(end+1) = size_Protein_data2; 
end
%
if length(data.se2) > 2
    data.se2 = data.se2(1:2);
end 

rng default % For reproducibility
data.Ansatz_          = rand(11,size_Protein_data+2);

end

