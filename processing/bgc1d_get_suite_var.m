function [varout] = bgc1d_get_suite_var(Suite,varname); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d ncycle v 6.1 - D. Bianchi 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get varaibale value for a suite, irrespective of suite and variable dimension
% size_out = [suite size]x[variable size]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if nargin==1
    varname = 'o2';
 end

 tmpvar = Suite.Out{1}.(varname);
 ssize = size(Suite.Out);	% suite size
 vsize = size(tmpvar);		% variable size

 nruns = prod(ssize);
 nval = prod(vsize);

 varout = nan(nruns,nval);

 % Reshapes Suite and variable so to create a 2D final output
 tmpOut = reshape(Suite.Out,1,nruns);

 % Creates indicing for Suite matrix
 for indr=1:nruns
    tmpvar = tmpOut{1,indr}.(varname);
    tmpvar = reshape(tmpvar,1,nval);
    varout(indr,:) = tmpvar; 
 end
 
 % Final reshape back to [ssize]x[vsize]
 varout = squeeze(reshape(varout,[ssize vsize]));  
 
 
