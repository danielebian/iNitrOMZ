function plotCMAESoptim(OptimPath,varargin)
A=[];
A=parse_pv_pairs(A,varargin)
% temporary hardcode to test function
OptimPath='/home/yangsi/project/NitrOMZ/iNitrOMZ_v6.0/comet/TestSY/optimize_cmaes/Optim_27_Nov_2019_03_22_06_FAOxAnox5_ni_cf3_nr_25k.mat';
rrpath='/home/yangsi/project/NitrOMZ/iNitrOMZ_v6.0/comet/TestSY/optimize_cmaes/';
addpath(rrpath);
addpath('/home/yangsi/project/NitrOMZ/optimization/CMA_ES/')

load(OptimPath);

% Structure is called Optim
load(OptimPath);
d.f=load([rrpath,'outcmaesfit.dat']);
d.x=load([rrpath,'outcmaesxmean.dat']);
d.sstd=load([rrpath,'outcmaesstddev.dat']);
d.bestever=load([rrpath,'outcmaesxrecentbest.dat']);


RunTimeHours=Optim.RunTime/3600;

rates={'KAo', 'KNo','KDen1','KDen2','KDen3','KAx'};
idx = zeros(size(Optim.ParNames));
for i = 1 : length(Optim.ParNames)
   if sum(strcmp(Optim.ParNames{i},rates))>0
      idx(i)=1;
   end
end
rateIdx=find(idx);
tracerIdx=find(~idx);

Optim.ParOptScaled=Optim.ParOpt;
Optim.ParOptScaled(rateIdx)=Optim.ParOpt(rateIdx)*86400;



%Plot Optimization performance
Fig=figure('units','centimeter','position',[0,0,17.9,12])
    ax(1)=subplot(2,2,1);
    ax(2)=subplot(2,2,2)
    ax(3)=subplot(2,2,3)
    ax(4)=subplot(2,2,4)
    subplot(ax(1))
        p(1)=semilogy(dfit.iter,dfit.mmin,'k','LineWidth',0.01);hold on;
        p(2)=semilogy(dfit.iter,dfit.mmax,'k','LineWidth',0.01);
        p(1).Color(4)=0.2;p(2).Color(4)=0.2;
        p(3)=semilogy(dfit.iter,dfit.best,'r');
        p(4)=semilogy(dfit.iter,dfit.sigma,'b');hold off;
        legend([p(3),p(4),p(1)],'Best','StdDev','min/max','Location','Northwest');
        ylabel('Cost')
        title(['Optim. convergence (',num2str(round(RunTimeHours),'%2d'),' hours)'])
    subplot(ax(2))
       plot(dparams.iter,dparams.mmean);
        ylim([0-0.2,1+0.2]);
        title('Parameters -- mean (normalized)');
    subplot(ax(3))
        bar(Optim.ParOptScaled(rateIdx));
        xticks((1:length(Optim.ParNames(rateIdx))));
        xticklabels(Optim.ParNames(rateIdx));
        xtickangle(45);
        title('Best params -- rates')
        ylabel('1/d')
    subplot(ax(4))
        bar(Optim.ParOptScaled(tracerIdx));
        xticks((1:length(Optim.ParNames(tracerIdx))));
        xticklabels(Optim.ParNames(tracerIdx));
        xtickangle(45);
        title('Best params -- tracers')
        ylabel('mmol/m3')
       
% Rerun model with Optimized parameters
