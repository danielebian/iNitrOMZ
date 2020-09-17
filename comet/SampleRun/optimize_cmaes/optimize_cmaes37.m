% Here, specify iNitrOMZ root path ($PATHTOINSTALL/iNitrOMZ/)
 bgc1d_root='/home/yangsi/project/NitrOMZ/iNitrOMZ_v6.0/';
% Adds CMAES subroutine:
 addpath('/home/yangsi/project/NitrOMZ/optimization/CMA_ES/');

% Handles output directory
 OptName = 'Optimmay437';
 OptCase = 'Optimmay4'
 OptNameLong = 'Fix All oxic, Vary Anoxic parameters, NO interpolated data, cost function #3 w/t no rates (&no second square), 25k evaluations max';
 Optim.codeDir = [bgc1d_root,'comet/',OptCase,'/optimize_cmaes/'];
 addpath(Optim.codeDir);
%OptName = 'Anox_Interp';
%OptName = 'Oxic_cmaes';
 curdir = pwd;
 DateNow = bgc1d_getDate();
% creates folder for CMAES output
 savedir = ['cmaes_out_' DateNow '_' OptName];
 mkdir([bgc1d_root 'optimOut'],savedir);
 cd([bgc1d_root 'optimOut/' savedir]);

% Parameters to tune:
 remin = 0.08/86400;

% Matrix of parameters for optimization
% Format: 	name 			min_value		max_value
 AllParam = {
               'wup_param',            0.4e-7,                 8.0e-7;         ...
%               'Kv_param',             0.5e-5,                 10e-5;          ...
%               'b',                    -1.0,                   -0.5;           ...
%                'poc_flux_top',         -15/86400,              -3/86400;       ...
%               'Krem',                 remin/10,               remin*5;        ...
                'Ji_a'                  0.05,                   0.4;          ...
                'Ji_b'                  0.05,                   0.2;          ...
                'KAo',                  0.05/86400,                    0.15/86400;              ...
                'KNo',                  0.05/86400,                  0.15/86400;              ...
                'KDen1',                remin/10,               remin;          ...
                'KDen2',                remin/10,               remin;          ...
                'KDen3',                remin/10,               remin;          ...
                'KAx',                  remin/5,               remin*5;        ...
                'KO2Rem',               0.01,                   1.0;            ...
                'KO2Den1',              0.01,                   6.0;            ...
                'KO2Den2',              0.01,                   3.0;            ...
                'KO2Den3',              0.01,                   3.0;              ...
                'KO2Ax',                0.5,                    6.0;              ...
                'KNH4Ao',               0.01,                   1.0;              ...
                'KNO2No',               0.01,                   1.0;            ...
                'KNO3Den1',             0.01,                   1.0;              ...
                'KNO2Den2',             0.01,                   1.0;            ...
                'KN2ODen3',             80/1000,               200/1000;            ...
                'KNH4Ax',               0.1,                   1.0;            ...
                'KNO2Ax',               0.1,                   1.0;            ...
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
 optn.EvalParallel = '1';
 optn.LBounds = (ParMin - ParMin) ./ ParNorm;
 optn.UBounds = (ParMax - ParMin) ./ ParNorm;
 optn.MaxFunEvals = 25000;

% Enables parallelization
% Note, the # of cores should be the same as the population size of the CMAES: 
% Popul size: 4 + floor(3*log(nPar))
 if strcmp(optn.EvalParallel,'1')
    FunName = 'bgc1d_fc2minimize_cmaes_parallel';
    delete(gcp('nocreate'))
    npar = 13;
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
 if strcmp(optn.EvalParallel,'1')
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
 Optim.bgc = bgc_run_Optim(Optim);

% Save ga output using today's date
 save(['Optim_' DateNow '_' OptName '.mat'],'Optim');
 cd(curdir)
