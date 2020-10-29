function Data = bgc1d_process_rates(Data1, bgc)

 % Initializes output
 Data = Data1;
 
 % Find index of o2 in solution "Sol"
 ind_o2 = find(strcmp(bgc.varname, 'o2'));
 
 % Case where Data only contains tracer data (i.e. first iteration);
 % Detects oxycline in the model, here defined at 1 mmol/m3
 depthox = bgc1d_detect_oxycline(bgc.sol(ind_o2,:),bgc);
 
 if isnan(depthox(1))
    % No oxycline. Assign large values for rates in order to kickout these runs.
    Data.rates.val = nan(length(Data.rates.name), length(bgc.zgrid));
    Data.rates.std = Data.rates.val;
    Data.rates.val(:) = 10^23;
    Data.rates.std(:) = 1.0;
 else
    Data.rates.val = nan(length(Data.rates.name), length(bgc.zgrid));
    Data.rates.std = Data.rates.val;
    % Here regrids observed rates onto model grid starting at the oxycline
    % Note that data rates are referenced to the oxycline using "depth_from_oxicline", z=0
    for indv=1:length(Data.rates.name)
       % Note conversion of observed rates from nM-N day^{-1} to mmol m^{-3} s^{-1}
       Data.rates.val(indv,:) = bgc1d_griddata(Data.rates.(Data.rates.name{indv}), ...
                                Data.rates.depth_from_oxicline+depthox(1),bgc) .*  ...
                                Data.rates.convf(indv); % convert from nM-N day^{-1} to mmol m^{-3} s^{-1}
    end
 end
 % Rearranges output structures
 if length(Data.name) < length(Data.name)+length(Data.rates.name)
    Data.val = vertcat(Data.val, Data.rates.val);
    Data.std = vertcat(Data.std, Data.rates.std);
    Data.weights = [Data.weights, Data.rates.weights];
 elseif length(Data.name)==length(Data.name)+length(Data.rates.name)
    Data.val(end-length(Data.rates.name):end,:) = Data.rates.val;
    Data.std(end-length(Data.rates.name):end,:) = Data.rates.std;
    Data.weights(end-length(Data.rates.name):end,:) = Data.rates.weights;
 else
    error('Data array is larger than combined Tracer and rate Array')
 end

