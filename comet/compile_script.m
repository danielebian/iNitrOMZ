%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to create compiled code for stand-alone matlab run
% Still in progress, but the basic idea is to create a working 
% directory locally where the code required is moved to, and
% and use it for compilation and running the new code.
% I clumsily address the problem of specifying function dependencies by
% making sure the main function paths are removed - the only path is local
% this is done by a separate function that is substituted before compilation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: 
% - Copy here the most up to date version sof the main optimization functions:
%   optimize_cmaes.m, bgc1d_fc2minimize_cmaes_parallel.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 curDir = pwd;

 rootPath = '/home/yangsi/project/NitrOMZ/';
 iNitPath = [rootPath,'iNitrOMZ_v6.0/'];
 addpath([iNitPath,'bgc1d_src/']);
 compDir = [iNitPath,'comet/']; 
 baseDir = 'Optimmay4';
 workDir = [compDir baseDir '/']; 

 % Directory for compilation
 mainFun = 'optimize_cmaes.m';
 % Function called by cmaes.m (cost function)
 auxFun = 'bgc1d_fc2minimize_cmaes_parallel.m';
 % Required data
 dataFile = 'compilation_ETSP_gridded_Feb232018_interpol.mat';
 funName = mainFun(1:end-2);
 compDir = [workDir funName];

 % Create directories for compilation
 if ~(exist(workDir)==7)
    mkdir(workDir)
 end
 if ~(exist(compDir)==7)
    mkdir(compDir)
 end

 % Main script to compile
 if (1)
    % Copy here the main functions required
    disp(['WARNING: copying main files locally']);
    copyfile([iNitPath,'runscripts/optimize_cmaes.m'],'.');
    copyfile([iNitPath,'optimization/bgc1d_fc2minimize_cmaes_parallel.m'],'.');
    copyfile([curDir,'/matlab.sb'],[workDir,'/optimize_cmaes/']);
    copyfile([curDir,'/runMultiMatlab.sb'],[workDir,'/optimize_cmaes/']);
 end
 sedcmd=['sed -i s/PpPlaceHoldeRrR/',baseDir,'/g ',workDir,'optimize_cmaes/','runMultiMatlab.sb'];
 system(sedcmd);
 mainFun = 'optimize_cmaes.m';
 % Function called by cmaes.m (cost function)
 auxFun = 'bgc1d_fc2minimize_cmaes_parallel.m';
 % Required data
 dataFile = 'compilation_ETSP_gridded_Feb232018_interpol.mat';

 % Compilation options
 % -m : specifies C language translation during compilation
 % -R : specifies runtime options
%compOpt = ['-m -R -singleCompThread -R -nosplash -R -nojvm']; 


 % Make sure local paths are present
 add_local_paths(rootPath)

 % Find file dependencies
 fList = matlab.codetools.requiredFilesAndProducts(mainFun);

 % Adds functions needed by auxFun
 fList = [fList matlab.codetools.requiredFilesAndProducts(auxFun)];

 % Adds any data needd
 fList = [fList which(dataFile)];

 fLisit = unique(fList);

 % Copies all file dependencies to compilation directory
 nList = length(fList);
 for indf=1:nList
   copyfile(fList{indf},compDir);
 end

 cd(compDir)

