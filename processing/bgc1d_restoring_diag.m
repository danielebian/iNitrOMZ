bgc1d_sms_diag.m function restoring = bgc1d_restoring_diag(bgc);

o2 = bgc.o2; no3 = bgc.no3; poc = bgc.poc; po4 = bgc.po4;
n2o=bgc.n2o; nh4 = bgc.nh4; no2 = bgc.no2; n2 = bgc.n2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d model v1.3 - Simon Yang  - October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies restoring source (linear relaxation) for selected dissolved tracers
 if bgc.tauZvar == 1
    tauh = ppval(bgc.restore.tauh_pre,bgc.zgrid);
 else
    tauh = bgc.tauh;
 end
 
 if bgc.PO4rest==0
    restoring.po4 = zeros(1,length(po4));
 else
    cout = ppval(bgc.restore.po4_cs,bgc.zgrid);
    restoring.po4 = (cout-po4)./tauh;
 end
 if bgc.NO3rest==0
    restoring.no3 = zeros(1,length(no3));
 else
    cout = ppval(bgc.restore.no3_cs,bgc.zgrid);
    restoring.no3 = (cout-no3)./tauh; 
 end
 if bgc.O2rest==0
    restoring.o2 = zeros(1,length(o2));
 else
    cout = ppval(bgc.restore.o2_cs,bgc.zgrid);
    restoring.o2 = (cout-o2)./tauh;
 end
 if bgc.N2Orest==0
    restoring.n2o = zeros(1,length(n2o));
 else
    cout = ppval(bgc.restore.n2o_cs,bgc.zgrid);
    restoring.n2o = (cout-n2o)./tauh;
 end
  if bgc.N2rest==0
    restoring.n2 = zeros(1,length(n2));
 else
    cout = zeros(size(bgc.zgrid));
    restoring.n2 = (cout-n2)./tauh;
  end
  if bgc.NO2rest==0
    restoring.no2 = zeros(1,length(no2));
 else
    cout = ppval(bgc.restore.no2_cs,bgc.zgrid);
    restoring.no2 = (cout-no2)./tauh;
  end
   if bgc.NH4rest==0
    restoring.nh4 = zeros(1,length(nh4));
 else
    cout = ppval(bgc.restore.nh4_cs,bgc.zgrid);
    restoring.nh4 = (cout-nh4)./tauh;
 end
