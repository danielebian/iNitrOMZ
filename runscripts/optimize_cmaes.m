 % Runs an optimization with CMA_ES approach
 bgc1d_root = '/Users/danielebianchi/AOS1/Ncycle/iNitrOMZ_v6.1/';
 %-----------------
 % Optimization names
 OptName = 'Opt_Nov_v1';
 OptNameLong = 'Fix baseline oxic param, vary Anoxic parameters, NO interpolated data, new rate constraints, *rescaled* by poc_flux ';

 imode=1; % 1 for interactive job; 2 for batch job

 switch imode
 case 1
    % Provides code paths:
    addpath([bgc1d_root 'bgc1d_src/']);
    addpath([bgc1d_root 'functions/']);
    addpath([bgc1d_root 'processing/']);
    addpath([bgc1d_root 'optimization/']);
    addpath([bgc1d_root 'runscripts/']);
    % Adds CMAES subroutine:
    addpath([bgc1d_root 'optimization/CMA_ES/']);
 case 2
   % When working on HPCS, move everything in one folder and specify path here:
   % Handles output directory using compilation script
   %addpath(Optim.codeDir);
 otherwise
    error(['Case not found']);
 end

 %-----------------------
 curdir = pwd;
% Creates folder for CMAES output
 DateNow = bgc1d_getDate();
 savedir = ['cmaes_out_' DateNow '_' OptName];
 mkdir([bgc1d_root 'optimOut'],savedir);
 cd([bgc1d_root 'optimOut/' savedir]);

%-------------------------------------------------------------------------------------
% Reference parameters
 remin = 0.08/86400;

% Parameters to tune:
% Matrix of parameters for optimization
% Format: 	name 			min_value		max_value
 AllParam = {
%               'wup_param',            0.4e-7,                 8.0e-7;         ...
%               'Kv_param',             0.5e-5,                 10e-5;          ...
%               'b',                    -1.0,                   -0.5;           ...
%               'poc_flux_top',         -15/86400,              -3/86400;       ...
%               'Krem',                 remin/10,               remin*5;        ...
                'Ji_a'                  0.05,                   0.4;            ...
                'Ji_b'                  0.05,                   0.2;            ...
                'KAo',                  0.01/86400,             0.50/86400;     ...
                'KNo',                  0.01/86400,             0.50/86400;     ...
                'KDen1',                remin/10,               remin;          ...
                'KDen2',                remin/10,               remin;          ...
                'KDen3',                remin/10,               remin;          ...
                'KAx',                  0.01/86400,             0.50/86400;     ...
                'KO2Rem',               0.01,                   1.0;            ...
                'KO2Den1',              0.01,                   6.0;            ...
                'KO2Den2',              0.01,                   3.0;            ...
                'KO2Den3',              0.01,                   3.0;            ...
                'KO2Ax',                0.5,                    6.0;            ...
                'KNH4Ao',               0.01,                   1.0;            ...
                'KNO2No',               0.01,                   1.0;            ...
                'KNO3Den1',             0.01,                   1.0;            ...
                'KNO2Den2',             0.01,                   1.0;            ...
                'KN2ODen3',             10/1000,               200/1000;        ...
                'KNH4Ax',               0.1,                   1.0;             ...
                'KNO2Ax',               0.1,                   1.0;             ...
                };

 ParNames = AllParam(:,1);
 nPar = size(ParNames,1); % number of parameters

 ParMin = [AllParam{:,2}]';
 ParMax = [AllParam{:,3}]';
 
 
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
%ParStart = (ParMean - ParMin) ./ ParNorm;
 ParStart = rand(nPar,1);
 ParSigma = ParRange./ParNorm/sqrt(12);

% Options
 optn.EvalParallel = 1;
 optn.LBounds = (ParMin - ParMin) ./ ParNorm;
 optn.UBounds = (ParMax - ParMin) ./ ParNorm;
 optn.MaxFunEvals = 25000;

% Enables parallelization
% Note, the # of cores should be the same as the population size of the CMAES: 
% Popul size: 4 + floor(3*log(nPar))
 if optn.EvalParallel==1
    FunName = 'bgc1d_fc2minimize_cmaes_parallel';
    delete(gcp('nocreate'))
    npar = 12;
    ThisPool = parpool('local',npar);
 else
    FunName = 'bgc1d_fc2minimize_cmaes';
 end

 FunArg.ParNames = ParNames;
 FunArg.ParMin = ParMin;
 FunArg.ParNorm = ParNorm;

% Runs the optimization
 tic;
 [pvarout, pmin, counteval, stopflag, out, bestever] = cmaes(FunName,ParStart,ParSigma,optn,FunArg);

% Stops parallel pool
 if optn.EvalParallel==1
    delete(ThisPool);
 end

% Fills in some output in final structure
% NOTE: instead of saving last iteration, saves best solution
 % Renormalized parameters
 Optim.OptNameLong = OptNameLong;
 Optim.ParOpt = ParMin + ParNorm .* bestever.x;
 Optim.ParNorm = ParNorm;
 Optim.cmaes.options = optn; 
 Optim.cmaes.pvarout = bestever.x; 
 Optim.cmaes.pmin = bestever.f; 
 Optim.cmaes.counteval = counteval; 
 Optim.cmaes.stopflag = stopflag; 
 Optim.cmaes.out = out; 
 Optim.cmaes.bestever = bestever; 
 Optim.RunTime = toc;
 % Runs and save best BGC1D parameters
 Optim.bgc = bgc1d_run('ParNames',Optim.ParNames,'ParVal',Optim.ParOpt);

% Save ga output using today's date
 save(['Optim_' DateNow '_' OptName '.mat'],'Optim');
 cd(curdir)
