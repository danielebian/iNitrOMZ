function bgc1d_plot(bgc,plot_mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d ncycle v 1.0 - Simon Yang  - October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results from optimization by the genetic algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot modelled O2, NO3, PO4, POC, N2O vs the Data used for the optimization%%%

 if nargin<2
    plot_mode = 1;
 end 

 % ------------------------------------------------
 % Plots O2, NO3, N2O
 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 10 7])

 %-------------
 subplot(2,2,1)
 s=scatter(bgc.Data_o2(~isnan(bgc.Data_o2)), bgc.zgrid(~isnan(bgc.Data_o2)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [0 0.5 0.5];
 grid on; hold on;
 plot(bgc.o2,bgc.zgrid,'k','linewidth',3)
 title(['O_2'] )
 ylabel('z (m)')
 xlabel('(mmol m^{-3})')
 ylim( [bgc.zbottom bgc.ztop]);
 xlim( [-5 max(bgc.o2)]);

 %-------------
 subplot(2,2,2)

 if (1)
    % Plots N* and N2 in N units
    plot(2*bgc.n2,bgc.zgrid,'--c','linewidth',3)
    grid on;hold on
    bgc.Data_nstar = bgc.NCrem/bgc.PCrem*bgc.Data_po4 - bgc.Data_no3;
    p=scatter(bgc.Data_nstar(~isnan(bgc.Data_nstar)), bgc.zgrid(~isnan(bgc.Data_nstar)),'b');
    p.MarkerEdgeColor = 'b';
    p.MarkerFaceColor = [0.1 0.1 0.4];
    plot(-bgc.nstar,bgc.zgrid,'-b','linewidth',3)
 end

 s=scatter(bgc.Data_no3(~isnan(bgc.Data_no3)), bgc.zgrid(~isnan(bgc.Data_no3)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [0 0.5 0.5];
 grid on;hold on
 plot(bgc.no3,bgc.zgrid,'-k','linewidth',3)
 plot(bgc.po4*bgc.NCrem/bgc.PCrem,bgc.zgrid,'--g','linewidth',3)
 
 p=scatter(bgc.Data_po4(~isnan(bgc.Data_po4))*bgc.NCrem/bgc.PCrem, bgc.zgrid(~isnan(bgc.Data_po4)),'b');
 p.LineWidth = 0.6;
 p.MarkerEdgeColor = 'g';
 p.MarkerFaceColor = [0 0.5 0.5];
 title(['NO_3, PO_4 ' ] )
 ylabel('z (m)')
 xlabel('(mmol m^{-3})')
 ylim( [bgc.zbottom bgc.ztop]);

 %-------------
 subplot(2,2,3)
 
 s=scatter(bgc.Data_nh4(~isnan(bgc.Data_nh4)), bgc.zgrid(~isnan(bgc.Data_nh4)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'r';
 s.MarkerFaceColor = [0 0.5 0.5];
 grid on; hold on;
 pp1 = plot(bgc.nh4,bgc.zgrid,'r','linewidth',3);
 ylabel('z (m)')
 xlabel('(mmol NH_4 m^{-3})')
 ylim( [bgc.zbottom bgc.ztop]);
 
 s=scatter(bgc.Data_no2(~isnan(bgc.Data_no2)), bgc.zgrid(~isnan(bgc.Data_no2)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [0 0.5 0.5];
 grid on; hold on;
 pp2 = plot(bgc.no2,bgc.zgrid,'k','linewidth',3);
 title(['NO_2, NH_4'] )
 ylabel('z (m)')
 xlabel('(mmol NO_2 m^{-3})')
 ylim( [bgc.zbottom bgc.ztop]);
 legend([pp1 pp2],'NH_4','NO_2','location','southwest');

 %-------------
 subplot(2,2,4)
 plot(bgc.n2o*1000,bgc.zgrid,'-k','linewidth',3)
 grid on; hold on
 s=scatter(bgc.Data_n2o(~isnan(bgc.Data_n2o))*1000, bgc.zgrid(~isnan(bgc.Data_n2o)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [0 0.5 0.5];
 title(['N_2O'] )
 ylabel('z (m)')
 xlabel('(nmol NO_2 m^{-3})')
 ylim( [bgc.zbottom bgc.ztop]);

 % More diagnostics:
 if plot_mode>1
 % ------------------------------------------------
 % Plots Remineralization diagnostics
 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 10 3.5])
 
 subplot(1,2,1)
 plot(bgc.RemOx./bgc.poc,bgc.zgrid,'k','linewidth',3)
 grid on;hold on;
 plot(bgc.RemDen./bgc.poc,bgc.zgrid,'r','linewidth',3)
 legend('Oxic','Denitrification')
 title(['Remineralization'] )
 ylabel('z (m)')
 xlabel('(s^{-1})')
 ylim( [bgc.zbottom bgc.ztop]);
 
 subplot(1,2,2)
 plot(bgc.RemDen./bgc.poc,bgc.zgrid,'k','linewidth',3)
 grid on; hold on;
 plot(bgc.RemDen1./bgc.poc,bgc.zgrid,'--b','linewidth',2)
 plot(bgc.RemDen2./bgc.poc,bgc.zgrid,'--g','linewidth',2)
 plot(bgc.RemDen3./bgc.poc,bgc.zgrid,'--r','linewidth',2)
 legend('Total','via NO_3 reduction','via NO_2 reduction','via N_2O reduction','Location','northeast')
 title(['Remin via Denitrification'] )
 ylabel('z (m)')
 xlabel('(s^{-1})')
 ylim( [bgc.zbottom bgc.ztop]);

 figure('units','inches','visible', bgc.visible)
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 10 7])
 
 % ------------------------------------------------
 % Plots Remineralization diagnostics
 subplot(2,2,1)
 plot(-bgc.NCden1*bgc.RemDen1+bgc.Nitrox,bgc.zgrid,'k','linewidth',3)
 grid on; hold on;
 plot(-bgc.NCden1*bgc.RemDen1,bgc.zgrid,'--b','linewidth',2)
 plot(bgc.Nitrox,bgc.zgrid,'--g','linewidth',2)
 legend('NO_3 SMS','NO_3 reduction','NO_2 oxidation','Location','southeast')
 title(['NO_3 sources and sinks'] )
 ylabel('z (m)')
 xlabel('(mmol N m^{-3} s^{-1})')
 ylim( [bgc.zbottom bgc.ztop]);

 subplot(2,2,2)
 plot(bgc.NCden1*bgc.RemDen1+bgc.Jno2_hx-bgc.NCden2*bgc.RemDen2-bgc.Nitrox-bgc.Anammox,bgc.zgrid,'k','linewidth',3)
 grid on; hold on;
 plot(-bgc.NCden2*bgc.RemDen2,bgc.zgrid,'--b','linewidth',2)
 plot(-bgc.Nitrox,bgc.zgrid,'--g','linewidth',2)
 plot(-bgc.Anammox,bgc.zgrid,'--r','linewidth',2)
 plot(bgc.Jno2_hx,bgc.zgrid,'--c','linewidth',2)
 plot(bgc.NCden1*bgc.RemDen1,bgc.zgrid,'--m','linewidth',2)
 legend('NO_2 SMS','NO_2 reduction','NO_2 oxidation','Anammox','NH_4 oxidation', ...
        'NO_3 reduction','Location','southeast')
 title(['NO_2 sources and sinks'] )
 ylabel('z (m)')
 xlabel('(mmol N m^{-3} s^{-1})')
 ylim( [bgc.zbottom bgc.ztop]);

 subplot(2,2,3)
 plot((0.5*bgc.NCden2*bgc.RemDen2-bgc.NCden3*bgc.RemDen3+bgc.Jnn2o_hx+bgc.Jnn2o_nden)*1000*3600*24,bgc.zgrid,'k','linewidth',3)
 grid on; hold on;
 plot((-bgc.NCden3*bgc.RemDen3)*1000*3600*24,bgc.zgrid,'--b','linewidth',2)
 plot((0.5*bgc.NCden2*bgc.RemDen2)*1000*3600*24,bgc.zgrid,'--g','linewidth',2)
 plot((bgc.Jnn2o_hx+bgc.Jnn2o_nden)*1000*3600*24,bgc.zgrid,'--r','linewidth',2)
 legend('N_2O SMS','N_2O reduction','NO_2 reduction','N_2O source from Ao','Location','southwest')
 title(['N_2O sources and sinks'] )
 ylabel('z (m)')
 xlabel('(nM N d^{-1})')
 ylim( [bgc.zbottom bgc.ztop]);
 
 subplot(2,2,4)
 plot(bgc.NCrem*(bgc.RemOx+bgc.RemDen1+bgc.RemDen2+bgc.RemDen3)-bgc.Ammox-bgc.Anammox,bgc.zgrid,'k','linewidth',3)
 grid on; hold on;
 plot(-bgc.Ammox,bgc.zgrid,'--b','linewidth',2)
 plot(-bgc.Anammox,bgc.zgrid,'--g','linewidth',2)
 plot(bgc.NCrem*(bgc.RemOx+bgc.RemDen1+bgc.RemDen2+bgc.RemDen3),bgc.zgrid,'--r','linewidth',2)
 legend('NH_4 SMS','NH_4 oxidation','Anammox','Reminerlization','Location','southwest')
 title(['NH_4 sources and sinks'] )
 ylabel('z (m)')
 xlabel('(mmol N m^{-3} s^{-1})')
 ylim( [bgc.zbottom bgc.ztop]);

 % ------------------------------------------------
 % Plots Flux diagnostics
 if bgc.flux_diag == 1
    figure('units','inches')
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 10 7])
    plot(bgc.advn2o, bgc.zgrid,'linewidth',3)
    grid on; hold on;
    plot(bgc.diffn2o, bgc.zgrid,'linewidth',3)
    plot(bgc.smsn2o, bgc.zgrid,'linewidth',3)
    plot(bgc.restn2o, bgc.zgrid,'linewidth',2)
    plot(bgc.advn2o+bgc.diffn2o+bgc.smsn2o+bgc.restn2o, bgc.zgrid,'linewidth',4)
    legend('N_2O advection','N_2O diffusion','SMS','N_2O restoring','Diffusion - advection + SMS + restoring')
    title(['N_2O advection, diffusion, SMS and restoring'] )
    ylabel('z (m)')
    xlabel('(mmol N2O m^{-2} s^{-1})')
    ylim( [bgc.zbottom bgc.ztop]);
 end
 
 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 10 7])
 % Use this for explicit sinking solution
%plot(bgc.poc.*bgc.wsink, bgc.zgrid,'--k','linewidth',3)
 % Use this for implicit sinking solution (flux shifted up by one level)
 plot([bgc.poc_flux_top bgc.poc(1:end-1).*bgc.wsink(1:end-1)], bgc.zgrid(1:end),'-k','linewidth',3)
 martin0p85 = @(z) bgc.poc_flux_top.*(z/bgc.ztop).^(-0.85);
 martin0p35 = @(z) bgc.poc_flux_top.*(z/bgc.ztop).^(-0.35);
 hold on;grid on;
 plot(martin0p85(bgc.zgrid), bgc.zgrid,'--r','linewidth',2)
 plot(martin0p35(bgc.zgrid), bgc.zgrid,'--g','linewidth',2)
 legend('POC flux','martin curve, b=-0.85','martin curve, b=-0.35')
 title(['N_2O advection, diffusion, SMS and restoring'] )
 ylabel('z (m)')
 xlabel('(mmol POC m^{-2} s^{-1})')
 end
 
 end

