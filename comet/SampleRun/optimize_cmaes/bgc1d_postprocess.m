 function bgc = bgc1d_postprocess(bgc, Data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d Ncycle v1.0 - Simon Yang  - April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the ini solution from the bgc.sol structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.sol = squeeze(bgc.sol_time(end,:,:));
 for indt=1:bgc.nvar
    bgc.(bgc.varname{indt}) = bgc.sol(indt,:);
    if bgc.flux_diag == 1
       bgc.(['adv' bgc.varname{indt}]) = bgc.sadv(indt,:);
       bgc.(['diff' bgc.varname{indt}]) = bgc.sdiff(indt,:);
       bgc.(['sms' bgc.varname{indt}]) = bgc.ssms(indt,:);
       bgc.(['rest' bgc.varname{indt}]) = bgc.srest(indt,:);
    end
    bgc.(['d' bgc.varname{indt}]) = nan(size(bgc.(bgc.varname{indt})));
    bgc.(['d' bgc.varname{indt}])(2:end-1) = ...
	(bgc.(bgc.varname{indt})(3:end) - bgc.(bgc.varname{indt})(1:end-2))/(-2*bgc.dz);
    bgc.(['d' bgc.varname{indt}])(1) = 0; 
    bgc.(['d' bgc.varname{indt}])(end) = 0;
    bgc.(['d2' bgc.varname{indt}]) = nan(size(bgc.(['d' bgc.varname{indt}])));
    bgc.(['d2' bgc.varname{indt}])(2:end-1) = ...
	(bgc.(bgc.varname{indt})(3:end) - 2 * bgc.(bgc.varname{indt})(2:end-1) ...
	+ bgc.(bgc.varname{indt})(1:end-2))/(bgc.dz^2);
    bgc.(['d2' bgc.varname{indt}])(1) = 0;
    bgc.(['d2' bgc.varname{indt}])(end) = 0;
 end
 if bgc.RunIsotopes
    ii = dNiso('i15N', bgc.i15no3, 'i14N', bgc.no3);
    idx = find(bgc.no3<bgc.IsoThreshold | bgc.i15no3<0);
    bgc.d15no3 = ii.d15N;
    %bgc.d15no3(idx)=nan;
    ii = dNiso('i15N', bgc.i15no2, 'i14N', bgc.no2);
    idx = find(bgc.no2<bgc.IsoThreshold | bgc.i15no2<0);
    bgc.d15no2 = ii.d15N;
    %bgc.d15no2(idx)=nan;
    ii = dNiso('i15N', bgc.i15nh4, 'i14N', bgc.nh4);
    idx = find(bgc.nh4<bgc.IsoThreshold | bgc.i15nh4<0);
    bgc.d15nh4 = ii.d15N;
    %bgc.d15nh4(idx)=nan;
    ii = dNiso('i15N', bgc.i15n2oA, 'i14N', bgc.n2o);
    idx = find(bgc.n2o<bgc.IsoThreshold/1000 | bgc.i15n2oA<0);
    bgc.d15n2oA = ii.d15N;
    %bgc.d15n2oA(idx)=nan;
    ii = dNiso('i15N', bgc.i15n2oB, 'i14N', bgc.n2o);
    idx = find(bgc.n2o<bgc.IsoThreshold/1000 | bgc.i15n2oB<0);
    bgc.d15n2oB = ii.d15N;
    %bgc.d15n2oB(idx)=nan;
 end

 if nargin>1
    ntrData = length(bgc.varname);
    tmp = strcat('Data_', bgc.varname);
    for indt=1:ntrData
       try
          bgc.(tmp{indt}) = Data.val(indt,:);
       catch
          display(['Warning: did not find', tmp{indt}]);
       end
    end
    try
       bgc.Data_nstar = bgc.Data_no3-bgc.NCrem/bgc.PCrem*bgc.Data_po4;
    end
 end

 % Additional tracers:
 % NSTAR - NO3 deficit versus PO4
 bgc.nstar = bgc.no3-bgc.NCrem/bgc.PCrem*bgc.po4;

 bgc.dwup = nan(size(bgc.wup)); bgc.dwup(1)=0; bgc.dwup(end)=0;
 bgc.dwup(2:end-1) = (bgc.wup(3:end) - bgc.wup(1:end-2))/(-2*bgc.dz);
 bgc.ssN2OAdv =  - (bgc.wup .* bgc.dn2o + bgc.dwup .* bgc.n2o);
 bgc.ssN2ODiff =bgc.Kv .* bgc.d2n2o; %bgc.vdKv .* bgc.dn2o + bgc.Kv .* bgc.d2n2o;
 bgc.sms_n2o(1)=0;bgc.sms_n2o(end)=0;bgc.n2o_rest(1) =0; bgc.n2o_rest(end)=0;
 % BGC sources and sink
 % Biological sources and sinks terms
 sms =  bgc1d_sms_diag(bgc); 
 bgc.RemOx    = sms.RemOx*1000*3600*24;
 bgc.Ammox    = sms.Ammox*1000*3600*24;
 bgc.Anammox  = sms.Anammox*1000*3600*24;
 bgc.Nitrox    = sms.Nitrox*1000*3600*24;
 bgc.RemDen   = sms.RemDen*1000*3600*24;
 bgc.RemDen1   = sms.RemDen1*1000*3600*24;
 bgc.RemDen2   = sms.RemDen2*1000*3600*24;
 bgc.RemDen3   = sms.RemDen3*1000*3600*24;
 bgc.Jnn2o_hx = sms.Jnn2o_hx*1000*3600*24; % nM N/d
 bgc.Jnn2o_nden = sms.Jnn2o_nden*1000*3600*24; % nM N/d
 bgc.Jno2_hx = sms.Jno2_hx*1000*3600*24; % nM N/d
 bgc.Jn2o_prod = sms.Jn2o_prod*1000*3600*24; % nM N2O/d
 bgc.Jn2o_cons = sms.Jn2o_cons*1000*3600*24;% nM N2O/d
 bgc.Jno2_prod = sms.Jno2_prod*1000*3600*24;% nM N/d
 bgc.Jno2_cons = sms.Jno2_cons*1000*3600*24;% nM N/d
 bgc.sms_n2o = sms.n2o*1000*3600*24; % nM N2O/d

 bgc.no2ton2o = sms.n2oind.den2;% nM N2O/d
 bgc.n2onetden = 0.5*bgc.NCden2*bgc.RemDen2-bgc.NCden3*bgc.RemDen3;% nM N2O/d
 bgc.nh4ton2o = 0.5*bgc.Jnn2o_hx+0.5*bgc.Jnn2o_nden;% nM N2O/d
 bgc.no3tono2 = bgc.NCden1*bgc.RemDen1;% nM N/d
 bgc.AnammoxFrac = (2*bgc.Anammox)./(2.*bgc.Anammox+bgc.NCden2.*bgc.RemDen2+bgc.Jn2o_prod.*2).*100;

 if bgc.RunIsotopes
    bgc.r15no3 = sms.r15no3;
    bgc.r15no2 = sms.r15no2;
    bgc.r15nh4 = sms.r15nh4;
    bgc.r15n2o = sms.r15n2o;
    bgc.r15n2oA = bgc.i15n2oA./bgc.n2o;
    bgc.r15n2oB = bgc.i15n2oB./bgc.n2o;
 end
