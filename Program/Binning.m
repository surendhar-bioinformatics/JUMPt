function data = binning(params)
%% read the data
All_data    = readtable(params.input_file);
time_points_ind 	= startsWith(All_data.Properties.VariableNames, 'time_point');
conc_ind            = startsWith(All_data.Properties.VariableNames, 'concentration');
data.t              = table2array(All_data(1,time_points_ind));  
data.t_long         = linspace(0, 32,321);
    
data.TotalLysConc	= table2array(All_data(2,conc_ind));
data.LysRatio    	= table2array(All_data(3,time_points_ind));% Lys data
data.freeLysConc    = table2array(All_data(3,conc_ind));
data.ProtInfo       = (All_data(4:end,1:find(time_points_ind,1)-1)); % protein IDs
data.ProtConc       = table2array(All_data(4:end,conc_ind))';% Protein concentration
SILAC_data          = table2cell(All_data(4:end,time_points_ind));% protein data
SILAC_data_NaN      = cellfun(@(x) ~isa(x,'double'),SILAC_data); % check for possible 0x0 char reads
SILAC_data(SILAC_data_NaN)   = {NaN};% set those to NaN
data.SILAC_data              = cell2mat(SILAC_data)'; % heavy/light ratios,
data.timepoints = All_data.Properties.VariableNames(time_points_ind);

% calculate apprent half-lives to sort the proteins
size_Prot_data  = size(data.SILAC_data,2);
fun_deg = @(b,tspan) exp(-b*tspan); 
for i =1:size_Prot_data 
    Pdata_temp  = data.SILAC_data(:,i);
	non_NaN_index_of_Pdata_temp = find(~isnan(Pdata_temp));
	Pdata_temp  = Pdata_temp(non_NaN_index_of_Pdata_temp);
    rng default % For reproducibility
	ExpDegRate_temp         = fminsearch(@(b) norm(Pdata_temp - fun_deg(b,data.t(non_NaN_index_of_Pdata_temp))), rand(1,1)) ;
	apparent_half_lives(i)  = ExpDegRate_temp;
end
apparent_half_lives = transpose(log(2)./apparent_half_lives);
[~, idx]            = sort(apparent_half_lives);
data.SILAC_data     = data.SILAC_data(:,idx); 
data.ProtConc       = data.ProtConc(:,idx);
data.ProtInfo       = data.ProtInfo(idx, :); 

% define Bin_size and Bin#
if size_Prot_data < params.bin_size
    params.bin_size = size_Prot_data;
end

% create bin array for all proteins and proteins with all time points 
if size_Prot_data == params.bin_size
    data.protBins           = [0 params.bin_size];
    data.protBins_allTimes  = [0 params.bin_size];%define bin for proteins with all time points
    data.SILAC_data_allTimes= data.SILAC_data;
    data.ProtInfo_allTimes  = data.ProtInfo;
    data.ProtConc_allTimes  = data.ProtConc;
else
    % create bin array for proteins with all time points 
    data.SILAC_data_allTimes= data.SILAC_data(:, (sum(~isnan(data.SILAC_data)) == sum(time_points_ind(:)== 1)));
    data.ProtInfo_allTimes  = data.ProtInfo((sum(~isnan(data.SILAC_data)) == sum(time_points_ind(:)== 1)),:);
    data.ProtConc_allTimes  = data.ProtConc(:,(sum(~isnan(data.SILAC_data)) == sum(time_points_ind(:)== 1)));
    i=1;
    while (size(data.SILAC_data_allTimes,2) < params.bin_size) && (i < sum(time_points_ind(:) == 1))
        data.SILAC_data_allTimes    = data.SILAC_data(:, (sum(~isnan(data.SILAC_data)) >= sum(time_points_ind(:) == 1)-i));
        data.ProtInfo_allTimes      = data.ProtInfo((sum(~isnan(data.SILAC_data)) >= sum(time_points_ind(:) == 1)-i),:);
        data.ProtConc_allTimes      = data.ProtConc(:, (sum(~isnan(data.SILAC_data)) >= sum(time_points_ind(:) == 1)-i));
        i=i+1;
    end
    bins        = ceil(size(data.SILAC_data_allTimes,2)/params.bin_size);
    indices     = ([0:params.bin_size-1]*bins)+1;
    indices_0   = indices;
    for i = 1:bins-1
        indices = [indices (indices_0+i)];
    end
    indices(indices > size(data.SILAC_data_allTimes,2)) = [];
    data.SILAC_data_allTimes= data.SILAC_data_allTimes(:,indices); 
    data.ProtConc_allTimes  = data.ProtConc_allTimes(:,indices);
    data.ProtInfo_allTimes  = data.ProtInfo_allTimes(indices, :);
    data.protBins_allTimes  = [0 params.bin_size];%creating an array with bin size

    % create bin array for all proteins 
    bins = ceil(size_Prot_data/params.bin_size);
    indices = ([0:params.bin_size-1]*bins)+1;
    indices_0 = indices;
    for i = 1:bins-1
        indices = [indices (indices_0+i)];
    end
    indices(indices > size_Prot_data) = [];
    data.SILAC_data = data.SILAC_data(:,indices); 
    data.ProtConc   = data.ProtConc(:,indices);
    data.ProtInfo   = data.ProtInfo(indices, :); 

    data.SILAC_data = data.SILAC_data(:, 1:size_Prot_data); 
    data.ProtConc   = data.ProtConc(:, 1:size_Prot_data);
    data.ProtInfo   = data.ProtInfo(1:size_Prot_data,:); 
    %creating an array with bin size
    data.protBins = [0];
    for i = 1:floor(size_Prot_data/params.bin_size)
        data.protBins(i+1) = i*params.bin_size;
    end
    if (mod(size_Prot_data,params.bin_size) > 0)  && (mod(size_Prot_data,params.bin_size) >= params.bin_size/4)
        data.protBins(end+1) = size_Prot_data; 
    elseif (mod(size_Prot_data,params.bin_size) > 0)  && (mod(size_Prot_data,params.bin_size) < params.bin_size/4) && (length(data.protBins) >1)
        data.protBins(end)   = size_Prot_data; 
    elseif (mod(size_Prot_data,params.bin_size) > 0)  && (mod(size_Prot_data,params.bin_size) < params.bin_size/4) && (length(data.protBins) == 1)
        data.protBins(end+1) = size_Prot_data; 
    end
end

rng default % For reproducibility
data.Ansatz_	= rand(11,size_Prot_data+2);% generate initial guess for the parametrs

end

