function addpath(rootpath,varargin)
% Opionally add more paths as a cell array {path1, path2, etc..}
 A.AddPaths =  [];
 A=parse_pv_pairs(A,varargin);
% Add paths for function dependencies
% Better be specific about the paths
 if ~isempty(A.AddPaths)
    for i = 1 : length(A.AddPaths)
       addpath(A.AddPaths{i});
    end
 end
 addpath([rootpath,'/optimization/CMA_ES/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/bgc1d_src/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/UserParams/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/functions/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/optimization/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/pprocessing/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/preprocessing/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/runscripts/']);
 addpath([rootpath,'/iNitrOMZ_v6.0/Data/']);

