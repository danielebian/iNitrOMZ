function bgc=getIsotopeDeltas(bgc)
% Calculates Isotope Deltas. The function assume that the tracers 
% have already been post-processed, i.e., the tracers are in the bgc 
% structure (bgc.o2, bgc.no3 etc...).

 ii = dNiso('i15N', bgc.i15no3, 'N', bgc.no3);
 idx = find(bgc.no3<bgc.IsoThreshold | bgc.i15no3<0);
 bgc.d15no3 = ii.d15N;
 bgc.d15no3(idx)=nan;
 ii = dNiso('i15N', bgc.i15no2, 'N', bgc.no2);
 idx = find(bgc.no2<bgc.IsoThreshold | bgc.i15no2<0);
 bgc.d15no2 = ii.d15N;
 bgc.d15no2(idx)=nan;
 ii = dNiso('i15N', bgc.i15nh4, 'N', bgc.nh4);
 idx = find(bgc.nh4<bgc.IsoThreshold | bgc.i15no2<0);
 bgc.d15nh4 = ii.d15N;
 bgc.d15nh4(idx)=nan;
 ii = dNiso('i15N', bgc.i15n2oA, 'N', bgc.n2o);
 idx = find(bgc.n2o<bgc.IsoThreshold/1000 | bgc.i15n2oA<0);
 bgc.d15n2oA = ii.d15N;
 bgc.d15n2oA(idx)=nan;
 ii = dNiso('i15N', bgc.i15n2oB, 'N', bgc.n2o);
 idx = find(bgc.n2o<bgc.IsoThreshold/1000 | bgc.i15n2oB<0);
 bgc.d15n2oB = ii.d15N;
 bgc.d15n2oB(idx)=nan;
