function Data1 = bgc1d_process_rates_opt(Data, bgc)

 Data1 = Data;

 % Fixes naming issue with older dataset, where "depth_from_oxicline" was used instead of "depth_from_oxycline"
 if ~isfield(Data1.rates,'depth_from_oxycline') & isfield(Data1.rates,'depth_from_oxicline')
    Data1.rates.depth_from_oxycline = Data1.rates.depth_from_oxicline;
 end

 % Find index of o2 in solution "Sol"
 ind_o2 = find(strcmp(bgc.varname,'o2'));

 % Case where Data only contains tracer data (i.e. first iteration, first chromosome)
 depthox = bgc1d_detect_oxycline(bgc.sol(ind_o2,:),bgc);

 if isnan(depthox(1))
    % No oxycline. Assign large values for rates in order to remove these runs.
    Data1.rates.val = nan(length(Data1.rates.name),length(bgc.zgrid));
    Data1.rates.val(:) = 10^23;
 else
    % Re-builds the observed rate profiles on the same grid of the model
    % Here uses the depth_from_oxicline to redistribute rates, using model's oxycline 
    % position as reference depth
    % Note rate conversion not strictly needed, all rates should be nM N/d (convf = 1)
    % However leaves the option to use different "convf" if needed
    Data1.rates.val = nan(length(Data1.rates.name),length(bgc.zgrid));
    for indv=1:length(Data1.rates.name)
       grid_data = bgc1d_griddata(Data1.rates.(Data1.rates.name{indv}),Data1.rates.depth_from_oxycline+depthox(1), bgc);
       Data1.rates.val(indv,:) = grid_data .* Data1.rates.convf(indv); % convert units if needed
    end
 end

 %---------------------------------------
 % In case of rescaling observed rates by a factor dependend on model/data POC flux
 % Calculates here this factor
 if Data.rescale_by_poc==1
    % Find depth index for reference flux depth
    indz = findin(abs(Data.poc_flux_ref_depth),abs(bgc.zgrid));
    % Adds estimate of particle flux
    % This is somewhat approximate because it's recalculated from POC
    % And sinking speed at the tracer cells
    poc_flux = -bgc.wsink .* bgc.poc*86400; % mmol C/m2/d
    poc_factor = abs(poc_flux(indz))/Data.poc_flux_ref; 
    Data1.rates.val = Data1.rates.val * poc_factor;       
 end
 %---------------------------------------

 if length(Data1.name) < (length(Data1.name) + length(Data1.rates.name))
    Data1.val = vertcat(Data1.val,Data1.rates.val);
    Data1.weights = [Data1.weights,Data1.rates.weights];
 elseif length(Data1.name) == (length(Data1.name) + length(Data1.rates.name))
    Data1.val(end-length(Data1.rates.name):end,:) = Data1.rates.val;
    Data1.weights(end-length(Data1.rates.name):end,:) = Data1.rates.weights;
 else
    error('Data array is larger than combined Tracer and rate Array');
 end

