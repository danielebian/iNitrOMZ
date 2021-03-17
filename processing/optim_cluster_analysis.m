function [idx creslt cdata] = optim_cluster_analysis(opt,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merges two optimizations, doing a KS test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
 A.nclust = 4;
 A.mode = 'kmeans';	% 'kmeans','gmm'
 A.iplot = 1;
 A.add_param = 1;
 A.add_cost = 0;
 A.add_var = 0;
 A.names_param = 'all';		% 'all' or 'none' or {'','',''}
%A.names_cost = [1 1 0 1 1 1 1 1 1 1 1 1 1];
 A.names_cost = [0 0 0 0 1 0 1 0 0 0 0 0 0];
 A.names_var = 'all';		% 'all' or 'none' or {'','',''}
 A.var_z_range = [40 600];
 A = parse_pv_pairs(A,varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Creates clustering dataset
 nrun = opt.nopt;

 cdata = nan(0,nrun); 

 % Adds parameters
 if A.add_param==1
    if ~iscell(A.names_param) & strcmp(lower(A.names_param),'all');
       parnames = opt.ParNames;
    elseif ~strcmp(lower(A.names_var),'none')
       parnames = A.names_param;
    else
       parnames = {};
    end

    % Create input dataset for clustering
    for indp=1:length(parnames)
       thisvar = opt.(parnames{indp}); 
       cdata = [cdata;thisvar'];
    end
 end

 % Adds parameters
 % varname: 'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2' 'nstar'
 % weights:  2    1     0     1     12     2     8     0    4
 % rates:   'nh4ton2o' 'noxton2o' 'no3tono2' 'anammox'
 % weights:  0          0          0          0
 if A.add_cost==1
    % First always adds the final cost. Here uses re-calculated 
    % on the original observational dataset (without randomness)
    cdata = [cdata;opt.cost_new'];
    % Then adds the variable-specific dataset, based on the vector "names_cost" of in (1) and out (0)
    % Note: order should be same as "bgc1d_fc2minimize_evaluate_cost.m", with variables and rates
    % Create input dataset for clustering
    tmpvar = opt.cost_pre(:,find(A.names_cost))';
    cdata = [cdata;tmpvar];
 end

 % Adds variables
 if A.add_var==1
    % Depth selection:
    zrange = -abs(A.var_z_range);
    % finds indices
    izrange = sort([findin(zrange(1),opt.var.zgrid) findin(zrange(2),opt.var.zgrid)]); 
    % Adds variables
    if ~iscell(A.names_param) & strcmp(lower(A.names_var),'all');
       varnames = fieldnames(opt.var);
       % Removes typically unwanted variables
       varnames = setdiff(varnames,{'zgrid','ncount','poc','AnammoxFrac','poc_flux'},'stable');
    elseif ~strcmp(lower(A.names_var),'none')
       varnames = A.names_var;
    else
       varnames = {};
    end
    for indv=1:length(varnames)
       thisvar = opt.var.(varnames{indv})(:,izrange(1):izrange(2))'; 
       cdata = [cdata;thisvar];
    end
 end

 % Performs kmeans clustering
 % First dataset needs to be of size NxP, with N number of instances, P of parameters for clustering
 cdata = cdata';

 switch A.mode
 case 'kmeans'
    [idx,C,sumd] = kmeans(cdata,A.nclust);
    rslt.C = C; 
    rslt.sumd = sumd; 
 case 'kmedoids'
    [idx,C,sumd] = kmedoids(cdata,A.nclust);
    rslt.C = C; 
    rslt.sumd = sumd; 
 case 'gmm'
    gm = fitgmdist(cdata,A.nclust);
    [idx,nlogL,P,logpdf,d2] = cluster(gm,cdata);
    rslt.nlogL = nlogL;
    rslt.P = P;
    rslt.nlogL = nlogL;
 case 'clusterdata'
    idx = clusterdata(cdata,'Maxclust',A.nclust);
 otherwise
    error(['Mode ' A.mode ' not supported']);
 end

 if A.iplot==1
    figure(987);
    [silh,h] = silhouette(cdata,idx);
    xlabel('Silhouette Value')
    ylabel('Cluster')
    title(['Mean silh: ' num2str(mean(silh))]);
 end








