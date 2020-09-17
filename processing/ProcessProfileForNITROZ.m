%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function grids  tracer profile data  for use with NITROMZ
%
% Simon Yang @ UCLA, June 6th 2019
%
% Assumes data is in .xlsx format, with each column being a different tracer (or depth/pressure etc...). 
%
% %%% Arguments:
% 	% xlsxPath ---> Path to excel file
%	% ColumnNames ---> Cell array containg the names of the variables corresponding
%			   to the columns in xlsxPath. e.g. {'Depth','O2' ,'NO3', 'PO4'}
% %%% Optional arguments:
%	% 'ztop' ---> Top of model grid in meters (see bgc1d_initialize)
%	% 'zbottom' ---> Bottom of model grid in meters (see bgc1d_initialize)
%	% 'npt' ---> Number of model grid points (see bgc1d_initialize)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


function compilation =  ProcessProfileForNITROZ(xlsxPath, ColumnNames, varargin)
% Default arguments
 A.ztop = -30; % top of model 
 A.zbottom = -1330;
 A.npt = 130;
 A = parse_pv_pairs(A,varargin);

% Calculate grid dimensions
 depth = abs(A.ztop-A.zbottom);               % Water column depth (m)
 nz=A.npt+1;                                  % Number of points in z
 dz = (A.ztop - A.zbottom)/A.npt;             % dz (m)
 zgrid = linspace(A.ztop,A.zbottom,A.npt+1);  % z points (m)
 zbounds = [zgrid+dz/2,zgrid(end)-dz/2]; % z bounds (m)

 vals = xlsread(xlsxPath) ;

 Depth=vals(:,1)*-1;
 for i = 2 : length(ColumnNames)-1
	 [smean] = gridprofile(Depth,vals(:,i), zgrid);
	 if strcmpi(ColumnNames{i},'n2o')
		 compilation.(ColumnNames{i}) = smean/1000;
	 else
	 	compilation.(ColumnNames{i}) = smean;
	end
 end

 compilation.zgrid=zgrid;

