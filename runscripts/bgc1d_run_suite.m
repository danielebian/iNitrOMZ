% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template BGC1D runscript 
% Versions: 0.1 : D. Bianchi, Spet 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Documentation:
% This script allows to perform "sensitivity experiments", where an arbitrary number
% of model parameters is varied, and the model run for each possible combinations of
% parameters.
% Output for each parameter combination will be stored in the cell array Suite.Out
% Access individual experiments by using the appropriate indices, e.g. Suite.Out{1,1}, ... 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Adds path for local functions

 %-------------------------------------------------------
 % Initialize the model - baseline
 % Based on the code and options in bgc1d_main.m
 %-------------------------------------------------------
 clear bgc;

 %-------------------------------------------------------
 bgc.root = '/Users/danielebianchi/AOS1/Ncycle/iNitrOMZ_v6.1/'; % root path
 bgc.rootdirs = {'bgc1d_src','Data','functions','processing','plotting','optimization','restart','runscripts'};
 bgc.iPlot = 0;            % To plot output
% initialize the model
% here, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
 for indp=1:length(bgc.rootdirs);
    addpath([bgc.root bgc.rootdirs{indp}]); % adds root to MATLAB's search paths
 end

 bgc.depparams = 1; % make sure dependent parameters are updated
 bgc = bgc1d_initialize(bgc);

 %-------------------------------------------------------
 % Process observations to validate the model solution
 Tracer.name = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'};
 if strcmp(bgc.region,'ETNP')
    load([bgc.root,'/Data/compilation_offshore.mat']);
    Data = GA_data_init_opt(bgc,compilation_offshore,Tracer.name);
 elseif strcmp(bgc.region,'ETSP')
   load([bgc.root,'/Data/compilation_ETSP_gridded_Feb232018.mat']);
   Data = GA_data_init_opt(bgc,compilation_ETSP_gridded,Tracer.name);
 end

 %-------------------------------------------------------
 % Here, if needed, overrides default parameters for BioModules and SetUp
 % All Suite experiments will adopt these parameters
 % (use ['property',value] format)
 % NOTE: these should be variables not used as Suite Parameters
%new_ParName = {'dt','nt','hist'};
%new_ParVal = [2*86400,250*365,365*10];
 new_ParName = {'dt','endTimey','histTimey'};
 new_ParVal = [2*86400,800,10];

 %-------------------------------------------------------
 % Define the suite of model runs
 %-------------------------------------------------------

 clear Suite;
 Suite.name = 'test';
 Suite.long_name = 'Test';
 NameAdd = 1;  % 1 to add the parameter names to the Suite name
 Suite.base = bgc;
 Suite.collapse = 1; 	% Collapses the suite Output by taking average output
			% and packaging the output into arrays with the size 
			% of the Suite parameters (useful to save space, removes time-dependent output)
 Suite.rmOut = 1;	% 0: keeps Out; 1: removes Out
 %---------------------
 % Suite parameters
 % Specify the following:
 % params : names of parameters to be varied
 % values : one vector of values for each parameter
 %---------------------
 Suite.params	= {'dt'};
%Suite.params	= {'poc_flux_top'};
%Suite.params	= {'wup_param','poc_flux_top'};

 Suite.dt = [2 4 8 16 32]*86400;
%Suite.wup_param = [0.25 0.5 1 2 4]*7.972e-8;
%Suite.poc_flux_top = [0.5:0.05:1.2]*(-7.5/86400); 

 %-------------------------------------------------------
 Suite.nparam = length(Suite.params);
 Suite.dims = zeros(1,Suite.nparam);
 Suite.AllParam = cell(1,Suite.nparam);
 for ip = 1:Suite.nparam
    Suite.dims(ip) = length(eval(['Suite.' Suite.params{ip}]));
    Suite.AllParam{ip} = eval(['Suite.' Suite.params{ip}]);
 end
 Suite.nruns = prod(Suite.dims);
 if length(Suite.dims)>1
    Suite.Out = cell(Suite.dims);
 else
    Suite.Out = cell(1,Suite.dims);
 end

 %-------------------------------------------------------
 % Loop through experiments
 %-------------------------------------------------------

 Tsuite = Suite.base;
 runindex = cell(Suite.nparam,1);
 for irun = 1:Suite.nruns
    disp(['Run number # ' num2str(irun) '/' num2str(Suite.nruns)]);
    [runindex{:}] = ind2sub(Suite.dims,irun);
    %---------------------
    % Here, creates the input arguments for bio module and experiment setup
    arg_ParName = new_ParName;
    arg_ParVal = new_ParVal;
    for ipar = 1:Suite.nparam
       disp([ Suite.params{ipar} ' - Start ........  ' num2str(Suite.AllParam{ipar}(runindex{ipar}))]);
       arg_ParName = [arg_ParName Suite.params{ipar}]; 
       arg_ParVal  = [arg_ParVal  Suite.AllParam{ipar}(runindex{ipar})]; 
    end
    %---------------------
    % Initialized Biomodule with user-defined inputs
    % (use ['property',value] format)
    % Substitute the optimal parameters
    % Change parameters with those selected in Scenario_child(ichr).chromosome
    for indp=1:length(arg_ParName)
       Tsuite = bgc1d_change_input(Tsuite,arg_ParName{indp},arg_ParVal(indp));
    end
    % Updates BGC/N-cycling parameters  that depend on bgc1d_initbgc_params
    if Tsuite.depparams
       Tsuite = bgc1d_initialize_DepParam(Tsuite);
       % Calculate dependent variables relate to isotopes
       if Tsuite.RunIsotopes
          Tsuite = bgc1d_initIso_Dep_params(Tsuite);
       end
    end
    % % % % % % % % % % % % % % % % % % % % 
    if (1)
       disp(['WARNING: updating any derived parameter']);
       % % % % % % % % % % % %
       % DERIVED PARAMETERS  %
       % % % % % % % % % % % %
       % E.g.  the following allows flexible timestepping
       Tsuite.nt = ceil(Tsuite.endTimey*(365*86400)/Tsuite.dt);
       Tsuite.hist = floor(Tsuite.histTimey*(365*86400)/Tsuite.dt);
       Tsuite.nhist = floor(Tsuite.nt/Tsuite.hist);
       Tsuite.hist_time = linspace(Tsuite.dt*Tsuite.hist, ...
                          Tsuite.nt*Tsuite.dt,Tsuite.nhist)/86400/365;
       % E.g. Constraints: 
       % ... 
       % ... KDen1 + KDen2 + KDen3 = remin
       % % % % % % % % % % % %
    end
    %---------------------
    Suite.Out{irun} = Tsuite;
    tic;
    % Run the model
    [Suite.Out{irun}.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(Suite.Out{irun});
    Suite.Out{irun}.runtime = toc;
    % Postprocess the results
    Suite.Out{irun} = bgc1d_postprocess(Suite.Out{irun}, Data);
    % Keeps track of runtime
 end
 %---------------------
 % Keeps track of total time, summing up individual times
 Suite.runtime = 0;
 for irun = 1:Suite.nruns
     Suite.runtime = Suite.runtime + Suite.Out{irun}.runtime;
 end;

 %-------------------------------------------------------
 % Postprocess, rename and save the suite
 %-------------------------------------------------------
 % If required, collapses Suite output
 if Suite.collapse==1
    % WARNING: this removes the "Out" field
 end

 % Rename the suite
 snewname = ['Suite_' Suite.name];
 if NameAdd ==1
    % Create a newname that includes all the parameters
    for indn=1:Suite.nparam
       snewname = [snewname '_' Suite.params{indn}];
    end
 end
 eval([snewname ' = Suite;']);
 % Save the suite
 eval(['save ' snewname ' ' snewname ';']);

