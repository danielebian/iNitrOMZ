function optout = optimization_add_cost(optin,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merges two optimizations, doing a KS test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
 A.iplot = 0;
 A = parse_pv_pairs(A,varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 optout = optin;
 optout.cost_new = nan(optout.nopt,1);
 % Hardwire # of cost constraints
 optout.cost_pre = nan(optout.nopt,13);

 for indr=1:optout.nopt
   tbgc = optout.bgc{indr};    
   [tcost cost_pre] = bgc1d_fc2minimize_evaluate_cost(tbgc,A.iplot);
   optout.cost_new(indr,1) = tcost;
   optout.cost_pre(indr,:) = cost_pre;
 end

