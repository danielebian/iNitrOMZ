function bgc1d_plot(bgc,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d ncycle v 1.0 - Simon Yang  - October 2017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results from optimization by the genetic algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: bgc1d_plot_var(bgc,'var',{'o2','po4','no3','no2'},'fig',3,'col',[1.0 0.5 0.5])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%A.var = 'o2';
 A.var = {'o2','po4','no3','no2','nh4','n2o','nstar','n2','poc'};
 A.fig = 0;
 A.col = [0 0 0];
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

 nvar = length(A.var);
 pp = numSubplots(nvar);
 sp1 = pp(1);
 sp2 = pp(2);

 for indv=1:nvar;

    varname = A.var{indv};
   
    var_plot = bgc.(varname); 
    if isfield(bgc,['Data_' varname])
       var_data = bgc.(['Data_' varname]);
       var_range = [min([var_plot,var_data]) max([var_plot,var_data])];
    else
       var_range = [min(var_plot) max(var_plot)];
    end
    var_span = diff(var_range);
   
    % ------------------------------------------------
    % Plots Variable
    subplot(sp1,sp2,indv)
    if isfield(bgc,['Data_' varname])
       s=scatter(var_data(~isnan(var_data)), bgc.zgrid(~isnan(var_data)),'b');
       s.LineWidth = 0.6;
       s.MarkerEdgeColor = 'k';
       s.MarkerFaceColor = [0.7 0.7 0.9];
       hold on; 
    end
    plot(var_plot,bgc.zgrid,'-','color',A.col,'linewidth',3)
    title([varname])
    ylabel('z (m)')
    xlabel([varname ' units'])
    ylim([bgc.zbottom bgc.ztop]);
    xlim([var_range(1)-var_span/10 var_range(2)+var_span/10]);
    grid on; box on; 
   
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
