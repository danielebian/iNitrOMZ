function bgc = bgc_run_Optim(Optim) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given an Optimal solution from the CMAES algorithm, 
% initializes and runs the BGC1D model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template iNitrOMZ runscript 
% Versions: 5.4 : Simon Yang  - April 2019
%           6.0 : Daniele Bianchi - June 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Customize your model run in bgc.root/UserParams/
%   % General model set-up 	 -- bgc1d_initialize.m
%   % Boundary conditions        -- bgc1d_initboundary.m
%   % BGC/N-cycling params       -- bgc1d_initbgc_params.m
%   % N-isotopes-cycling params  -- bgc1d_initIso_params.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ParNames = Optim.ParNames;
 ParVal = Optim.ParOpt;

% HERE, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
 bgc1d_root='/u/home/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/';
 addpath(genpath(Optim.codeDir)); % adds root to MATLAB's search paths

% initialize the model
 clear bgc;
 bgc.depparams = 1; % make sure Dependant parameters are updated
 bgc = bgc1d_initialize(bgc); 

% Substitute the optimal parameters
% Change parameters with those selected in Scenario_child(ichr).chromosome
 for indp=1:length(ParNames)
    bgc = bgc1d_change_input(bgc,ParNames{indp},ParVal(indp));
 end

 if (1)
    % % % % % % % % % % % %
    % For consistency, NEED to follow what is done in:
    % bgc1d_fc2minimize_cmaes_parallel.m
    % % % % % % % % % % % %
    % DERIVED PARAMETERS  %
    % % % % % % % % % % % %
    % Here imposes the following constraint:
    % Constraints: KDen1 + KDen2 + KDen3 = remin
   %remin = 0.08/86400;
   %bgc.KDen1 = max(0, remin - bgc.KDen2 - bgc.KDen3);
   %disp(['WARNING: overriding "KDen1=remin-bgc.KDen2-bgc.KDen3" to match optimization constraints']);
    % % % % % % % % % % % %
 end


 % Update dependent parameters
 bgc = bgc1d_initialize_DepParam(bgc);
 if bgc.RunIsotopes
    bgc = bgc1d_initIso_Dep_params(bgc);
 end

% run the model 
% % % % % % % % % 
%     % bgc.sol_time is a TxVxZ matrix where T is archive times
%     	  V is the number of tracers and Z in the number of 
%	  model vertical levels
%     % note that the model saves in order:
%	  (1) o2 (2) no3 (3) poc (4) po4 (5) n2o (6) nh4 (7) no2 (8) n2
% 	  (9) i15no3 (10) i15no2 (11) i15nh4 (12) i15n2oA (13) i15n2oB 
 tic;
 [bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
 bgc.sol = squeeze(bgc.sol_time(end,:,:));
 bgc.RunTime = toc;
 disp(['Runtime : ' num2str(bgc.RunTime)]);

% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % Below are optional routines % % % % % 
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
 rates.name = {'nh4ton2o' 'noxton2o' 'no3tono2' 'anammox'} 
 [~,~,bgc] = bgc1d_fc2minimize_evaluate_cost(bgc);
 nconst=size(bgc.constraints_data,1)
 for i = 1 : length(rates.name)
    bgc.(['Data_',rates.name{i}]) = bgc.constraints_data(nconst-length(rates.name)+1,:);
 end
% Process model output for analysis (gathers tracers and diagnostics into the bgc structure)
 bgc = bgc1d_postprocess(bgc, Data);
 if (0)
    bgc1d_plot(bgc); 
 end

