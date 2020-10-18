function bgc = bgc1d_getrates(bgc,Data)

 zstep = bgc.npt;
 bgc.z = linspace(bgc.zbottom,bgc.ztop,zstep+1);
 bgc.dz = (bgc.ztop - bgc.zbottom) / zstep;
 bgc.SolNames = {'o2','no3','poc','po4','n2o', 'nh4', 'no2', 'n2'};
 ntrSol = length(bgc.SolNames);
 for indt=1:ntrSol
    bgc.(bgc.SolNames{indt}) = bgc.sol(indt,:);
 end

%error(['THIS NEEDS TO BE UPDATED TO USE SOURCESINK INSTEAD OF SMS_DIAG']);
%[sms diag] =  bgc1d_sms_diag(bgc);

 % dump tracers in a structure "tr" - one by one (avoids eval)
 tr.o2  = bgc.o2;
 tr.no3 = bgc.no3;
 tr.poc = bgc.poc;
 tr.po4 = bgc.po4;
 tr.n2o = bgc.n2o;
 tr.nh4 = bgc.nh4;
 tr.no2 = bgc.no2;
 tr.n2  = bgc.n2;
 if bgc.RunIsotopes
    tr.i15no3  = bgc.i15no3;
    tr.i15no2  = bgc.i15no2;
    tr.i15nh4  = bgc.i15nh4;
    tr.i15n2oA = bgc.i15n2oA;
    tr.i15n2oB = bgc.i15n2oB;
 end

 [sms diag] =  bgc1d_sourcesink(bgc,tr);

 bgc.rates = nan(length(Data.rates.name),length(bgc.zgrid));
 tmp.nh4tono2 = diag.Ammox;
 tmp.anammox  = diag.Anammox;
 tmp.nh4ton2o = sms.n2oind.ammox + sms.n2oind.nden;
 tmp.noxton2o = bgc.NCden2*diag.RemDen2;
 tmp.no3tono2 = bgc.NCden1*diag.RemDen1;
 for indr = 1:length(Data.rates.name)
    bgc.rates(indr,:) = tmp.(Data.rates.name{indr});
 end
