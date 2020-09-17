function Data1 = bgc1d_process_rates_opt(Data, bgc)

 Data1 = Data;

 % Find index of o2 in solution "Sol"
 ind_o2 = find(strcmp(bgc.varname,'o2'));

 % Case where Data only contains tracer data (i.e. first iteration, first chromosome)
 depthox = bgc1d_detect_oxycline(bgc.sol(ind_o2,:),bgc);

 if isnan(depthox(1))
    % No oxycline. Assign large values for rates in order to kickout these runs.
    Data1.rates.val = nan(length(Data1.rates.name),length(bgc.zgrid));
    Data1.rates.val(:) = 10^23;
 else
    Data1.rates.val = nan(length(Data1.rates.name),length(bgc.zgrid));
    for v=1:length(Data1.rates.name)
       Data1.rates.val(v,:) = bgc1d_griddata(Data1.rates.(Data1.rates.name{v}),Data1.rates.depth_from_oxicline+depthox(1), bgc) .* ...
                              Data1.rates.convf(v); % convert from nM-N day^{-1} to mmol m^{-3} s^{-1}
    end
 end

 if length(Data1.name) < (length(Data1.name) + length(Data1.rates.name))
    Data1.val = vertcat(Data1.val,Data1.rates.val);
    Data1.weights = [Data1.weights,Data1.rates.weights];
 elseif length(Data1.name) == (length(Data1.name) + length(Data1.rates.name))
    Data1.val(end-length(Data1.rates.name):end,:) = Data1.rates.val;
    Data1.weights(end-length(Data1.rates.name):end,:) = Data1.rates.weights;
 else
    error('Data array is larger than combined Tracer and rate Array');
 end



			

