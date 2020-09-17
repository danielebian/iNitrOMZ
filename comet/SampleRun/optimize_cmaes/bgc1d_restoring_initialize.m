 function  bgc = bgc_initialize_restoring(bgc);
 %------------------------------------------------
 % Interpolate profile by creating a cubic spline 
 % NOTE: modify substituting [0;tconc;0] if flat-slope endings are needed
 % WARNING: avoid extrapolation by using profiles from 0-4000 m (even if sparse)
 %------------------------------------------------

 load([bgc.root,bgc.farfield_profiles]);
 load([bgc.root,bgc.Tau_profiles]);
 %------------------------------------------------
 bgc.restore = [];
 if  strcmp(bgc.region, 'ETNP')
    Farfield = farfield_ETNP_gridded;
 elseif strcmp(bgc.region, 'ETSP')
    Farfield = farfield_ETSP_gridded;
 end
 
 if  bgc.tauZvar == 1
    rest.tauh = ((rest.currentZ_etnp./(bgc.Lh) + rest.kappaZ_etnp.*2./((bgc.Lh)^2)).^-1);
    tdepth = -abs(rest.depth);
    tconc  = rest.tauh;     
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];     
    bgc.restore.tauh_pre = spline(tdepth,tconc);
 end

 if ~(bgc.PO4rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.po4;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.po4_cs = spline(tdepth,tconc);
    bgc.restore.po4_cout = ppval(bgc.restore.po4_cs,bgc.zgrid);    
 end

 if ~(bgc.NO3rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.no3;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.no3_cs = spline(tdepth,tconc);
    bgc.restore.no3_cout = ppval(bgc.restore.no3_cs,bgc.zgrid);    
 end

 if ~(bgc.O2rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.o2;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.o2_cs = spline(tdepth,tconc);
    bgc.restore.o2_cout = ppval(bgc.restore.o2_cs,bgc.zgrid);    
 end

 if ~(bgc.N2Orest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.n2o;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.n2o_cs = spline(tdepth,tconc);
    bgc.restore.n2o_cout = ppval(bgc.restore.n2o_cs,bgc.zgrid);    
 end

 if ~(bgc.NO2rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.no2;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.no2_cs = spline(tdepth,tconc);
    bgc.restore.no2_cout = ppval(bgc.restore.no2_cs,bgc.zgrid);    
 end

 if ~(bgc.N2rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.n2;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.n2_cs = spline(tdepth,tconc);
    bgc.restore.n2_cout = ppval(bgc.restore.n2_cs,bgc.zgrid);
 end

 if ~(bgc.i15NO3rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.i15no3;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.i15no3_cs = spline(tdepth,tconc);
    bgc.restore.i15no3_cout = ppval(bgc.restore.i15no3_cs,bgc.zgrid);
 end

 if ~(bgc.i15NO2rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.i15no2;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.i15no2_cs = spline(tdepth,tconc);
    bgc.restore.i15no2_cout = ppval(bgc.restore.i15no2_cs,bgc.zgrid);
 end

 if ~(bgc.i15NH4rest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.i15nh4;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.i15nh4_cs = spline(tdepth,tconc);
    bgc.restore.i15nh4_cout = ppval(bgc.restore.i15nh4_cs,bgc.zgrid);
 end

 if ~(bgc.i15N2OArest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.i15n2oA;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.i15n2oA_cs = spline(tdepth,tconc);
    bgc.restore.i15n2oA_cout = ppval(bgc.restore.i15n2oA_cs,bgc.zgrid);
 end 

 if ~(bgc.i15N2OBrest==0)
    tdepth = -abs(Farfield.zgrid);
    tconc  = Farfield.i15n2oB;
    ibad = find(isnan(tconc));
    tdepth(ibad) = [];
    tconc(ibad) = [];
    bgc.restore.i15n2oB_cs = spline(tdepth,tconc);
    bgc.restore.i15n2oB_cout = ppval(bgc.restore.i15n2oB_cs,bgc.zgrid);
 end

