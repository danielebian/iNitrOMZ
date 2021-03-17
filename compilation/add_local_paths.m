function add_local_paths(rootpath,varargin)
% Opionally add more paths as a cell array {path1, path2, etc..}
 if nargin==1
    AddPaths =  {};
 else
    AddPaths = varargin{1};
 end
% Add paths for function dependencies
% Better be specific about the paths
 if ~isempty(AddPaths)
    for i = 1 : length(AddPaths)
       addpath(AddPaths{i});
    end
 end
 addpath([rootpath,'bgc1d_src/']);
 addpath([rootpath,'functions/']);
 addpath([rootpath,'optimization/']);
 addpath([rootpath,'optimization/CMA_ES/']);
 addpath([rootpath,'processing/']);
 addpath([rootpath,'runscripts/']);
 addpath([rootpath,'restart/']);
 addpath([rootpath,'Data/']);

