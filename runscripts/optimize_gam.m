% Here, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
 bgc1d_root='/u/home/d/danieleb/NitrOMZ/iNitrOMZ_v6.0/';
 addpath(genpath(bgc1d_root)); % adds root to MATLAB's search paths
% Adds CMAES subroutine:
 addpath('/u/home/d/danieleb/NitrOMZ/optimization/CMA_ES');

% Handles output directory
%OptName = 'Anox_Interp';
 OptName = 'Anoxic_cf2';
 curdir = pwd;
 DateNow = bgc1d_getDate();
% creates folder for CMAES output
 savedir = ['cmaes_out_' DateNow '_' OptName];
 mkdir([bgc1d_root 'optimOut'],savedir);
 cd([bgc1d_root 'optimOut/' savedir]);

% Parameters to tune:
 remin = 0.08/86400;

% Matrix of parameters for optimization
% Format:       name                    min_value               max_value
 AllParam = {
%               'wup_param',            0.4e-7,                 8.0e-7;         ...
%               'Kv_param',             0.5e-5,                 10e-5;          ...
                'poc_flux_top',         -15/86400,              -3/86400;       ...
                'Krem',                 remin/10,               remin*5;        ...
                'KO2Rem',               0.01,                   5.0;            ...
%               'b',                    -1.0,                   -0.5;           ...
                'KNO2No',               0.01,                   0.5;            ...
                'KO2Den1',              0.01,                   10.0;            ...
                'KO2Den2',              0.01,                   10.0;            ...
                'KO2Den3',              0.01,                   10.0;              ...
                'KNO3Den1',             0.01,                   10.0;              ...
%               'KDen1',                remin/10,               remin*5;          ...
                'KDen2',                remin/10,               remin*5;          ...
                'KDen3',                remin/10,               remin*5;          ...
                'KNO2Den2',             0.001,                  0.5;            ...
                'KN2ODen3',             0.001,                  0.5;            ...
                'KAx',                  remin/10,               remin*5;        ...
                'KNH4Ax',               0.01,                   0.5;            ...
                'KNO2Ax',               0.01,                   0.5;            ...
                'KO2Ax',                0.01,                   10.0;              ...
                'KAo',                  remin/10,               remin*5;              ...
                'KNH4Ao',               0.001,                  2;              ...
                'KO2Ao',                0.001,                  10.0;              ...
                };

 ParNames = AllParam(:,1);
 ParMin = [AllParam{:,2}]';
 ParMax = [AllParam{:,3}]';

 nPar = size(ParNames,1); % number of parameters

% Initialize final output structure
 Optim.ParNames = ParNames;
 Optim.ParMin = ParMin;
 Optim.ParMax = ParMax;

% NOTES: 
% (1) Parameters are normalized by subtracting the min and dividing by
%     a normalization factor, typically the range (so they are b/w 0-1)
%     This is done to allow the CMAES algorithm to work in the space [0 1]
% (2) If needed, remember to add the constraint:
%     Constraints: KDen1 + KDen2 + KDen3 = remin
%     This should be done in the cost function (bgc1d setup step)
%     as an ad-hoc constraint (removes one degree of freedom)
%     Remember to remove the corresponding K from the parameter pool!
%     (suggestion: remove KDen1, since first step drives everuthing)
% Calculates useful quantities for normalization, optimization, etc.
 ParMean = (ParMin + ParMax)/2';
 ParRange = ParMax - ParMin;
 ParNorm = ParRange;
 ParStart = (ParMean - ParMin) ./ ParNorm;
 ParSigma = ParRange./ParNorm/sqrt(12);

% Options
 optn.EvalParallel = 1;
 optn.LBounds = (ParMin - ParMin) ./ ParNorm;
 optn.UBounds = (ParMax - ParMin) ./ ParNorm;

% Constraints: KDen1 + KDen2 + KDen3 = remin
% NOTE: This constraint can be included in the GAM automatically
% But has to be hard-coded cor the CMAES optimization algorithm
%Aeq = [0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0];
%Aeq = [    0 0 0 0 0 0 1 1 1 0 0 0 0 0 0];
%beq = remin;
 
 FunArg.ParNames = ParNames;
 FunArg.ParMin = ParMin;
 FunArg.ParNorm = ParNorm;

% Makes handle for cost function. x is an array of parameter values, param are the names which need to be passed 
 costfunc = @(x)bgc1d_fc2minimize_gam(x,FunArg);

% Options
 options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf, ...
                        'UseParallel', optn.EvalParallel, 'UseVectorized', false,'Display','iter');

 if optn.EvalParallel 
    delete(gcp('nocreate'))
    npar = 7;
    ThisPool = parpool(npar);
 end

% Runs the optimization
 tic
 param_out = ga(costfunc,nPar,[],[],[],[],optn.LBounds,optn.UBounds,[],options);

% Stops parallel pool
 if strcmp(optn.EvalParallel,'1')
    delete(ThisPool);
 end

% Fills in some output in final structure
% NOTE: instead of saving last iteration, saves best solution
 % Renormalized parameters
 Optim.ParOpt = ParMin + ParNorm .* param_out;
 Optim.ParNorm = ParNorm;
 Optim.ga.options = options;
 Optim.ga.options_setup = optn;
 Optim.RunTime = toc;
 % Runs and save best BGC1D parameters
 Optim.bgc = bgc_run_Optim(Optim);

% Save ga output using today's date
 save(['Optim_' DateNow '_' OptName '.mat'],'Optim');
 cd(curdir)
