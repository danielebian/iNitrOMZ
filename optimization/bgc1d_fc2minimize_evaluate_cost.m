function [tcost cost_pre bgc] = bgc1d_fc2minimize_evaluate_cost(bgc,iplot)

 if nargin<2
    iplot = 0;
 end

 % Evaluate cost function for individual bgc_realization

 % Tracer data
 % RAW DATA:
 tmp = load([bgc.root,'/Data/compilation_ETSP_gridded_Feb232018.mat']);
 % INTERPOLATED (SMOOTH) DATA:
%tmp = load([bgc.root,'/Data/compilation_ETSP_gridded_Feb232018_interpol.mat']);

 Data.name = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2','nstar'};
 Data = GA_data_init_opt(bgc,tmp.compilation_ETSP_gridded,Data.name);

 % IF NSTAR is a tracer used for the comparison, then adds it to bgc:
 if any(strcmp(Data.name,'nstar'))
    % tracers = {'o2', 'no3','poc', 'po4', 'n2o', 'nh4', 'no2', 'n2'};
    bgc.tracers = [bgc.tracers 'nstar'];
    bgc.sol = [bgc.sol ; 16 * bgc.sol(4,:) - bgc.sol(2,:) ];
 end

 % Rate data (will be processed online to adjust for oxycline depth)
 tmp = load([bgc.root,'/Data/comprates_ETSP.mat']);
 Data.rates = tmp.comprates.ETSP;
 Data.rates.name = {'nh4ton2o' 'noxton2o', 'no3tono2', 'anammox'};
 Data.rates.convf = [1/1000/(3600*24), 1/1000/(3600*24), 1/1000/(3600*24) , 1/1000/3600];
 
 %specify weights of tracer data for the optimization
%        bgc.varname = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2' 'nstar'}
%Data.weights =  	[2    0     0     1     0     0     0     0    0];	% Oxic optimization
%Data.weights =  	[1    1     0     1     2     0     1     0    0];	% Anoxic Optimization-1
 Data.weights =  	[2    0     0     1     6     0     3.0   0    2];	% Anoxic Optimization-1
 %			'nh4ton2o' 'noxton2o' 'no3tono2' 'anammox'
%Data.rates.weights = 	[1          1          1          1];
%Data.rates.weights = 	[0          0          0          0];
 Data.rates.weights = 	[1          1          0          1];
 
 if sum(Data.rates.weights) > 0
    bgc = bgc1d_getrates(bgc, Data);
    constraints_model = vertcat(bgc.sol, bgc.rates);
    constraints_data = bgc1d_process_rates_opt(Data, bgc);
    % Set rate data in top 50 m (3 cells) to 0
    depthNoRates = 50;
    idx = find(-bgc.zgrid<=depthNoRates);
    constraints_data.val(length(Data.weights)+1:end,idx) = nan;
 else
    constraints_model = bgc.sol;
    constraints_data = Data;
 end

 bgc.constraints_model = constraints_model;
 bgc.constraints_data = constraints_data.val;
 % If needed, uses depth-dependent weights
 % The vertical profiles will be weighted by the selected weights
 % 1. No depth-dependent weigths (all depths weighted the same, weights=1)
 % 2. Gaussian centered at OMZ core weighted up to 50% more
 % 3. Double gaussian -- upper and lower weighted up to 50% more
 idepth_weights = 2;
 switch idepth_weights
 case 1
    depth_weights = ones(1,size(constraints_model,2));
 case 2
    % Weights a gaussian center on OMZ 50% more
    zlev = [1:size(constraints_model,2)];
    zmeanloc = 20;
    zmeanstd = 15;
    depth_weights = 0.5 + 0.5*exp(-(zlev-zmeanloc).^2/(2*zmeanstd^2));
 case 3 
    % Weights a gaussian center on OMZ 50% more
    zlev = [1:size(constraints_model,2)];
    zmeanloc1 = 12;
    zmeanloc2 = 35;
    zmeanstd1 = 4;
    zmeanstd2 = 10;
    depth_weights = 0.5 + 0.5*exp(-(zlev-zmeanloc1).^2/(2*zmeanstd1^2)) +...
                        + 0.5*exp(-(zlev-zmeanloc2).^2/(2*zmeanstd2^2));
 otherwise
    error(['Need a valid option for depth weigths']);
 end

% *********************************************************************
% -- Evaluation --
% *********************************************************************
% Calculate Cost
% Note - should weight NaNs in bgc1d output heavily, to avoid 
% the corresponding parameter combinations
 cost_pre = Cost_quad_feb22_v2(constraints_model,constraints_data,depth_weights,iplot);

 % calculate mean of costs
 cost_pre(isnan(cost_pre)) = 0;
 % Need to add costs quadratically to prevent any one constraints to drift to far off.
 tcost = nansum((cost_pre.* constraints_data.weights').^2)/(nansum(Data.weights)+nansum(Data.rates.weights));
 %tcost = nansum(cost_pre.* constraints_data.weights')/(nansum(Data.weights)+nansum(Data.rates.weights));
 
