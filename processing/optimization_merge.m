function opt3 = optimization_merge(opt1,opt2,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merges two optimizations, doing a KS test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
 A.placeholder = [];    % placeholder variable
 A = parse_pv_pairs(A,varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 varadd = {'cost','counteval','runtime'}; 
 [varnames ii1 ii2] = intersect(opt1.ParNames,opt2.ParNames,'stable'); 
 
 opt3.indir = '';
 opt3.suffix = [opt1.suffix '_' opt2.suffix];
 opt3.nopt = opt1.nopt + opt2.nopt;
 opt3.names = [opt1.names;opt2.names];
 opt3.ParNames = varnames;
 npar = length(opt3.ParNames);
 opt3.ParMin = nan(npar,1);
 opt3.ParMax = nan(npar,1);
 opt3.OptNameLong = ['Merged ' opt1.suffix ' and ' opt2.suffix];
 for indv=1:npar
    mmin = min(opt1.range.(opt3.ParNames{indv})(1),opt2.range.(opt3.ParNames{indv})(1));
    mmax = max(opt1.range.(opt3.ParNames{indv})(2),opt2.range.(opt3.ParNames{indv})(2));
    opt3.range.(opt3.ParNames{indv}) = [mmin mmax];
 end
 opt3.ParOpt = [opt1.ParOpt(:,ii1);opt2.ParOpt(:,ii2)];
 opt3.ParNorm = [opt1.ParNorm(:,ii1);opt2.ParNorm(:,ii2)];
 allnames = [varadd(:);opt3.ParNames];
 for indv=1:length(allnames)
    opt3.(allnames{indv}) = [opt1.(allnames{indv});opt2.(allnames{indv})]; 
 end
 opt3.bgc = [opt1.bgc;opt2.bgc];
 % Here assumes optimizations have the same vertical grid
 opt3.var.zgrid = opt1.var.zgrid;
 opt3.var.ncount = [1:opt3.nopt];
 vname1 = fieldnames(opt1.var);
 vname2 = fieldnames(opt2.var);
 vname3 = intersect(vname1,vname1,'stable');
 vname3 = setdiff(vname3,{'zgrid','ncount'},'stable');
 for indv=1:length(vname3)
    opt3.var.(vname3{indv}) = [opt1.var.(vname3{indv});opt2.var.(vname3{indv})];
 end

 % Performs a KS test to check that the parameters come from the same distribution
 % 0: same distribution
 % 1: different distribution
 opt3.kstest = nan(npar,1);
 for indv=1:npar
   tks = kstest2(opt1.(varnames{indv}),opt2.(varnames{indv}));
   opt3.kstest(indv) = tks;
 end


