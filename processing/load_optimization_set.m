function optOut = load_optimization_set(varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads results from a series of optimizations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 
% optOut = load_optimization_set('indir','/Users/danielebianchi/AOS1/Ncycle/iNitrOMZ_v6.1/optimOut/from_hoffman2/optimOut/','suffix','Opt_Nov_v3_nir40k');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
 A.indir = '/Users/danielebianchi/AOS1/Ncycle/iNitrOMZ_v6.1/optimOut/from_hoffman2/optimOut/';
 A.suffix = 'Opt_Nov_v3_nir40k';
 A.varprof = {'o2','no3','poc','po4','n2o','nh4','no2','n2','nstar', ...
              'remox','ammox','anammox','nitrox','remden1','remden2', ...
              'remden3','nh4tono2','no2tono3','no3tono2','no2ton2o', ...
              'n2oton2','noxton2o','n2onetden','nh4ton2o','AnammoxFrac','poc_flux'};
 A = parse_pv_pairs(A,varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 indir = A.indir;
 suffix = A.suffix;

% Looks at contents of indir folder
 if ~strcmp(indir(end),'/')
    indir = [indir '/'];
 end

 % Assumes all folders end with "suffix"
 dirData = dir([indir '*' suffix]);
 if isempty(dirData)
    error(['Folder ' indir '*' suffix ' not found']);
 end

 dirnames = {dirData.name};
 dirnames = dirnames(:);

 ndirs = length(dirnames);

 % Initialize optimization set structure
 rawOut = struct;
 allNames = {};

 % Loops through all files
 for indf=1:ndirs
    disp(['Processing folder # ' num2str(indf) '/' num2str(ndirs)]);
    % Looks for an optimization result "Optim_*" within the folder
    tmpName = dir([indir dirnames{indf} '/Optim_*']);
    if ~isempty(tmpName)
      tmpFile = load([tmpName.folder '/' tmpName.name]); 
      tOptim = tmpFile.Optim;
      rawOut.Optim{indf} = tOptim;
      allNames{indf} = tmpName.name; 
    else
       disp(['No Optim file found in folder... ' dirnames{indf} '/']);
    end
 end
 % Removes empty structures from rawOut
 ibad = cellfun(@isempty,rawOut.Optim); 
 rawOut.Optim(ibad) = [];
 allNames(ibad) = [];
 nOut = length(rawOut.Optim);

 % Processes Optimization structures
 clear optOut
 optOut.indir = indir;
 optOut.suffix = suffix;
 optOut.nopt = nOut;
 optOut.names = allNames(:);
 optOut.ParNames = rawOut.Optim{1}.ParNames;
 optOut.ParMin = rawOut.Optim{1}.ParMin;
 optOut.ParMax = rawOut.Optim{1}.ParMax;
 optOut.OptNameLong = rawOut.Optim{1}.OptNameLong;
 nvar = length(optOut.ParNames);
 % Adds in variable min and max as range
 for indv=1:nvar
    optOut.range.(optOut.ParNames{indv}) = [optOut.ParMin(indv,1) optOut.ParMax(indv,1)];
 end
 % Initializes variables
 optOut.ParOpt = nan(nOut,nvar);
 optOut.ParNorm = nan(nOut,nvar);
 optOut.cost = nan(nOut,1);
 optOut.counteval = nan(nOut,1);
 optOut.runtime = nan(nOut,1);
 for indv=1:nvar
    optOut.(optOut.ParNames{indv}) = nan(nOut,1);
 end
 % Fills in all variables 
 for indo=1:nOut
    optOut.ParOpt(indo,:) = rawOut.Optim{indo}.ParOpt';
    optOut.ParNorm(indo,:) = rawOut.Optim{indo}.ParNorm';
    optOut.cost(indo) = rawOut.Optim{indo}.cmaes.pmin;
    optOut.counteval(indo) = rawOut.Optim{indo}.cmaes.counteval;
    optOut.runtime(indo) = rawOut.Optim{indo}.RunTime;
    optOut.bgc{indo,1} = rawOut.Optim{indo}.bgc;
    for indv=1:nvar
       optOut.(optOut.ParNames{indv})(indo) = rawOut.Optim{indo}.bgc.(optOut.ParNames{indv});
    end 
 end
 
 % If needed, creates 2D structures of final profiles from a list of provided
 % variables that are contained in the output structure "Optim.bgc"
 if ~isempty(A.varprof)
    nvp = length(A.varprof);
    nz = optOut.bgc{1,1}.nz;
    % Initializes
    optOut.var.zgrid = optOut.bgc{1,1}.zgrid;
    optOut.var.ncount = [1:nOut]; 
    for indv=1:nvp
       optOut.var.(A.varprof{indv}) = nan(nOut,nz);
    end
    % Fills in profiles for selected variables
    for indo=1:nOut
       for indv=1:nvp
          optOut.var.(A.varprof{indv})(indo,:) = optOut.bgc{indo,1}.(A.varprof{indv});
       end
    end
 end




















  



