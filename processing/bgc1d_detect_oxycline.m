  function depth = bgc1d_detect_oxycline(o2,bgc);

 depth = nan(1,2);

 ind_anoxix = find(o2 < 1.0);

 if length(ind_anoxix) >= 2
    depth(1) = bgc.zgrid(ind_anoxix(1));
    depth(2) = bgc.zgrid(ind_anoxix(end));
 elseif length(ind_anoxix) == 1
    depth(1) = bgc.zgrid(ind_anoxix(1));
    depth(2) = nan;
 elseif length(ind_anoxix) == 0
    depth(1) = nan;
    depth(2) = nan;
 end

