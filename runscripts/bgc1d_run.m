function bgc = bgc1d_run(varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template iNitrOMZ runscript 
% Versions: 5.4 : Simon Yang  - April 2019
%           6.0 : Simon Yang - June 2020
%           6.1 : Daniele Bianchi - September 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Customize your model run in bgc.root/UserParams/
%   % General model set-up 	 -- bgc1d_initialize.m
%   % Boundary conditions        -- bgc1d_initboundary.m
%   % BGC/N-cycling params       -- bgc1d_initbgc_params.m
%   % N-isotopes-cycling params  -- bgc1d_initIso_params.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A.iPlot = 0;		% To plot output
A.ParNames = {};	% Pass parameter names that need to be modified from default values
A.ParVal = [];		% Pass parameter values, corresponding to ParNames
A.derived_inputs = 0;	% Allows specific changes to variables by deriving selected parameters
			% from other parameters. This may be required for example for optimizations
A = parse_pv_pairs(A, varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the model
 clear bgc;
 bgc = struct;
% here, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
%bgc.root = A.root;
%for indp=1:length(A.rootdirs);
%   addpath([bgc.root A.rootdirs{indp}]); % adds root to MATLAB's search paths
%end

%bgc.depparams = 1; % make sure dependent parameters are updated
 bgc = bgc1d_initialize(bgc); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% In case parameters are specified as inputs, e.g. by passing the
% results of an Optimization then substitutes parameters
% ParNames = {'par1';'par2';...}
% ParVal = [val1; val2; ...];
% % % % % % % % % % % % % % % % % % % % 
% Substitute the optimal parameters
% Change parameters with those selected in Scenario_child(ichr).chromosome
 if ~isempty(A.ParNames) 
    for indp=1:length(A.ParNames)
       bgc = bgc1d_change_input(bgc,A.ParNames{indp},A.ParVal(indp));
    end
    % Updates BGC/N-cycling parameters  that depend on bgc1d_initbgc_params
    if bgc.depparams
       bgc = bgc1d_initialize_DepParam(bgc);
       % Calculate dependent variables relate to isotopes
       if bgc.RunIsotopes
          bgc = bgc1d_initIso_Dep_params(bgc);
       end
    end
 end
% % % % % % % % % % % % % % % % % % % % 
 if (A.derived_inputs)
    % % % % % % % % % % % %
    % For consistency, NEED to follow what is done in:
    % bgc1d_fc2minimize_cmaes_parallel.m
    % % % % % % % % % % % %
    % DERIVED PARAMETERS  %
    % % % % % % % % % % % %
    % Here imposes any required constraint:
    % E.g. Constraints: KDen1 + KDen2 + KDen3 = remin
    % remin = 0.08/86400;
    % bgc.KDen1 = max(0, remin - bgc.KDen2 - bgc.KDen3);
    % disp(['WARNING: overriding "KDen1=remin-bgc.KDen2-bgc.KDen3" to match optimization constraints']);
    % % % % % % % % % % % %
 end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % % % % % % % % % % % % % % % % % % % % 
% run the model 
% % % % % % % % % % % % % % % % % % % % % % % % 
%     % bgc.sol_time is a TxVxZ matrix where T is archive times
%     	  V is the number of tracers and Z in the number of 
%	  model vertical levels
%     % note that the model saves in order:
%	  (1) o2 (2) no3 (3) poc (4) po4 (5) n2o (6) nh4 (7) no2 (8) n2
% 	  (9) i15no3 (10) i15no2 (11) i15nh4 (12) i15n2oA (13) i15n2oB 
% % % % % % % % % % % % % % % % % % % % % % % % 
 tic;
 [bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
 bgc.RunTime = toc;
 disp(['Runtime : ' num2str(bgc.RunTime)]);

% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % Below are optional routines % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % %


% % % % % % % % % % % % % % % % % % % % % % % % 
% Alternatively, the model can output both tracers and all their fluxes 
% % % % % % % % %
%     % User must specify bgc.flux_diag == 1; in initialization function
%     % bgc.sadv_time + bgc.sdiff_time + bgc.ssms_time 
%          + bgc.srest_time =d(bgc.sol_time)/dt
%
%if bgc.flux_diag == 1 
%	[bgc.sol_time bgc.adv_time bgc.diff_time bgc.sms_time bgc.rest_time] = bgc1d_advection_diff(bgc);
%	bgc.sol = squeeze(bgc.sol_time(end,:,:)); % solution
%	bgc.sadv = squeeze(bgc.adv_time(end,:,:)); % advective fluxes
%	bgc.sdiff = squeeze(bgc.diff_time(end,:,:)); % diffusive fluxes
%	bgc.ssms = squeeze(bgc.sms_time(end,:,:)); % sources minus sinks
%	bgc.srest = squeeze(bgc.rest_time(end,:,:));  % restoring fluxes
%end

% Process observations to validate the model solution
 Tracer.name = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'};
 if strcmp(bgc.region,'ETNP')
    load([bgc.root,'/Data/compilation_offshore.mat']);
    Data = GA_data_init_opt(bgc,compilation_offshore,Tracer.name);
 elseif strcmp(bgc.region,'ETSP')
   load([bgc.root,'/Data/compilation_ETSP_gridded_Feb232018.mat']);
   Data = GA_data_init_opt(bgc,compilation_ETSP_gridded,Tracer.name);
 end

% Process model output for analysis (gathers tracers and diagnostics into the bgc structure)
 bgc = bgc1d_postprocess(bgc, Data);
 if (A.iPlot)
    bgc1d_plot(bgc); 
 end

