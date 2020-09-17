function bgc = bgc1d_getTracers(bgc)
        for indt=1:bgc.nvar
                bgc.(bgc.varname{indt}) = bgc.sol(indt,:);
        end
