function optim_plot_histograms(opt,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot histograms for parameters from an optimization
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
 A.var = {};
 A.num = 20;	% # bins for histograms
 A.mode = 1; 	% 1: linear space plot, prescribed # of bins
		% 2: log space, prescribed # of bins
		% 3: automatic histogram binning
 A.Normalization = 'count';	% help histogram, use 'count','probability', etc.
 A.BinMethod = 'auto';	% help histogram
 A.fig = 0;
 A = parse_pv_pairs(A,varargin);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if A.fig==0
   hfig = figure;
 else
   hfig = figure(A.fig);
 end

 if ~iscell(A.var)
    A.var = {A.var};
 end
 
 if isempty(A.var)
    A.var = opt.ParNames;
 end

 nvar = length(A.var);
 pp = numSubplots(nvar);
 sp1 = pp(1);
 sp2 = pp(2);

 for indv=1:nvar

    varname = A.var{indv};
    tvar = opt.(varname);
    varrange = opt.range.(varname);

    sp = subplot(sp1,sp2,indv);
    switch A.mode
    case 1
       % Linear range
       val0 = varrange(1);
       val1 = varrange(2); 
       hEdges = linspace(val0,val1,A.num);      
       hh = histogram(tvar,hEdges,'Normalization',A.Normalization);
    case 2
       val0 = log10(varrange(1));
       val1 = log10(varrange(2)); 
       hEdges = linspace(val0,val1,A.num);      
       hh = histogram(log10(tvar),hEdges,'Normalization',A.Normalization);
       xtk = get(sp,'XTick');
       xtk1 = 10.^xtk;
       xtk1 = mat2cell(xtk1,1,ones(size(xtk1)))';
       xlb = cellfun(@(x) num2str(x,2),xtk1,'UniformOutput',0);
       set(sp,'XTick',xtk,'XTickLabel',xlb);
    case 3
       hh = histogram(tvar,'Normalization',A.Normalization,'BinMethod',A.BinMethod);
    otherwise
       error(['Mode not found']);
    end 
    ysp = get(sp,'YLim');
    set(sp,'YLim',[ysp(1) ysp(2)*1.1]);

    xlabel(varname);

 end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots. 
%
% Inputs
% n - the desired number of subplots.     
%  
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.     
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p = 
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
%
%
% Rob Campbell - January 2010


while isprime(n) & n>4,
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design 
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
