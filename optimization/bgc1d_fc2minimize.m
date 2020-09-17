function cost = bgc1d_fc2minimize(parval,parnames)

 % Original list of parameters tested by Simon:
 % {'KO2Rem','KNO2No' 'KO2Den1','KO2Den2' 'KO2Den3','KNO3Den1','KDen1','KDen2', ...
 %  'KDen3','KNO2Den2','KN2ODen3','KAx','KNH4Ax','KNO2Ax','KO2Ax'};

 % Initialize model:
 bgc = bgc1d_initialize;
 bgc.hist_verbose = false;

 % Change parameters with those selected in Scenario_child(ichr).chromosome
 for indp=1:length(parnames)
    bgc = bgc1d_change_input(bgc,parnames{1,indp},parval(indp));
 end

 % update dependent parameters
 bgc = bgc1d_initialize_DepParam(bgc);
 if bgc.RunIsotopes
    bgc = bgc1d_initIso_Dep_params(bgc);
 end

 % % % % % % %
 % Run model %
 % % % % % % %
 [bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
 bgc.sol = squeeze(bgc.sol_time(end,:,:));
 
 % Tracer data
 load([bgc.root,'/Data/compilation_ETSP_gridded_Feb232018.mat']);
 Data.names = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'};
 Data = GA_data_init_opt(bgc,eval(['compilation_ETSP_gridded']),Data.names);

 % Rate data (will be processed online to adjust for oxycline depth)
 load([bgc.root,'/Data/comprates_ETSP.mat']);
 Data.rates = comprates.ETSP;
 Data.rates.name = {'nh4ton2o' 'noxton2o', 'no3tono2', 'anammox'};
 Data.rates.convf = [1/1000/(3600*24), 1/1000/(3600*24), 1/1000/(3600*24) , 1/1000/3600];
 
 %specify weights of tracer data for the optimization
 %'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'
 Data.weights =  [1  1  0  1  2  0  1  0];
 %'nh4ton2o' 'noxton2o' 'no3tono2' 'anammox'
 Data.rates.weights = [1  1  1  1];
 
 if sum(Data.rates.weights) > 0
    bgc = bgc1d_getrates(bgc, Data);
    data_constraints = vertcat(bgc.sol, bgc.rates);
    Data_processed = bgc1d_process_rates_opt(Data, bgc);
 else
    Data_processed = Data;
    data_constraints = bgc.sol;
 end

% *********************************************************************
% -- Evaluation --
% *********************************************************************
%   Calculate Cost
 cost_pre = Cost_quad_feb22_v2(data_constraints, Data_processed);

 % calculate mean of costs
 cost_pre(isnan(cost_pre)) = 0;
 cost = nansum((cost_pre.* Data_processed.weights').^2)/(nansum(Data.weights)+nansum(Data.rates.weights));

