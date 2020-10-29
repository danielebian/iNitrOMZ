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

 % NEW:
 % Converts from (uM N/s) to (nM N/d)
 % For safety recalculates and keeps rate into temporary structure "tmp"
 cnvrt = 1000*3600*24;
 tmp.nh4tono2   =       diag.Ammox   * cnvrt;      	% nM n/d
 tmp.anammox    = 2.0 * diag.Anammox * cnvrt;   	% nM N/d : Units of N, not N2
 tmp.no2tono3   =       diag.Nitrox  * cnvrt;      	% nM n/d
 tmp.nh4ton2o   = (diag.Jnn2o_hx + diag.Jnn2o_nden) * cnvrt; 	% nM N/d : Units of N, not N2O
 tmp.no3tono2   = bgc.NCden1 * diag.RemDen1 * cnvrt;	% nM n/d
 tmp.no2ton2o   =   2 * sms.n2oind.den2 * cnvrt;	% nM N/d : Units of N, not N2O
 tmp.n2oton2    = - 2 * sms.n2oind.den3 * cnvrt;	% nM N/d : Units of N, not N2
 tmp.noxton2o   = tmp.no2ton2o;				% nM N/d : Units of N, not N2O

 % Folds rates inside structure for easy handling in optimization
 bgc.rates = nan(length(Data.rates.name),length(bgc.zgrid));
 for indr = 1:length(Data.rates.name)
    bgc.rates(indr,:) = tmp.(Data.rates.name{indr});
 end

 % OLD:
 % Note, from "bgc1d_sourcesink.m" rates are returned in mmol/m3/s
%tmp.nh4tono2 = diag.Ammox;					% mmolN/m3/s 
%tmp.anammox  = diag.Anammox					% mmolN2/m3/s
%tmp.nh4ton2o = sms.n2oind.ammox + sms.n2oind.nden; 		% mmolN2O/m3/s
%tmp.noxton2o = bgc.NCden2*diag.RemDen2;			% mmolN/m3/s
%tmp.no3tono2 = bgc.NCden1*diag.RemDen1;			% mmolN/m3/s
