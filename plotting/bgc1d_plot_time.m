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
 A.var = {'o2','po4','no3','no2','nh4','n2o','poc','n2'};
 A.mode = 1;	% 1 plots variable
		% 2 plots anomaly versus final value
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

 plot_time = bgc.hist_time;

 for indv=1:nvar;

    varname = A.var{indv};
    isol = find(strcmp(varname,bgc.varname)); 

    if ~isempty(isol)
       switch A.mode
       case 1
          % plots variable
          var_plot = squeeze(bgc.sol_time(:,isol,:))';
       case 2
          % plots variable - end value
          var_plot = squeeze(bgc.sol_time(:,isol,:))';
          var_plot = var_plot - var_plot(:,end);
       case 3
          % plots variable - end value
          var_plot = squeeze(bgc.sol_time(:,isol,:))';
          var_plot = var_plot - var_plot(:,1);
       case 4
          % plots variable - end value
          var_plot = squeeze(bgc.sol_time(:,isol,:))';
          var_plot = var_plot - nanmean(var_plot,2);
       otherwise
         error('Plotting mode not found');
       end
       %var_range = [min(var_plot(:)) max(var_plot(:))];
       eps = 1e-10;
       mmin = min(var_plot(:));
       mmax = max(var_plot(:));
       var_range = [(1-sign(mmin)*0.1)*mmin-eps (1+sign(mmax)*0.1)*mmax+eps];
       var_span = diff(var_range);

       % ------------------------------------------------
       % Plots Variable
       subplot(sp1,sp2,indv)
       sanePColor(bgc.hist_time,bgc.zgrid,var_plot);
       shading flat
       colorbar
       title([varname])
       ylabel('depth (m)')
       xlabel('time (y)')
       ylim([bgc.zbottom bgc.ztop]);
       if A.mode==1
          caxis([var_range(1) var_range(2)]);
       else
          absmax = max(abs(mmin),abs(mmax));
          caxis([-absmax-eps absmax+eps]);
          colormap(mycolormaps('num',30,'scheme','redblue'));
       end
       grid on; box on; 
    end
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
