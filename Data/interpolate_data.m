
 if (1)
    addpath(genpath('/data/project1/matlabpathfiles/SLMtools/'));
    addpath('/data/project1/matlabpathfiles/');
    tmp = load('compilation_ETSP_gridded_Feb232018.mat')
    comp = tmp.compilation_ETSP_gridded;
    clear tmp
    comp.nstar = 16*comp.po4 - comp.no3;
 end

 zgrid = comp.zgrid;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % O2: Split into 3 regions:
 % From inspection of O2 and NO2:
 odz_up = 110;
 odz_lo = 360;
 iodz_up = findin(odz_up,zgrid);
 iodz_lo = findin(odz_lo,zgrid)-3;

 comp.o2_int = nan(size(comp.o2)); 

 % (1) oxic upper:
 tmpz = zgrid(1:iodz_up);
 tmpv = comp.o2(1:iodz_up);

 slm_presc = slmset('degree',3, ...
             'knots', 3, ...
             'minvalue',0, ...
             'decreasing','on', ...
             'concaveup', 'on', ...
             'rightvalue', 0, ...
             'rightslope', 0, ...
             'plot','off' ...
             ); 
 slm = slmengine(tmpz,tmpv,slm_presc);

 comp.o2_int(1:iodz_up) = slmeval(tmpz,slm);

 % (2) anoxic :
 comp.o2_int(iodz_up+1:iodz_lo-1) = 0;

 % (3) oxic lower:
 tmpz = zgrid(iodz_lo:end);
 tmpv = comp.o2(iodz_lo:end);

 slm_presc = slmset('degree',3, ...
             'minvalue',0, ... 
             'knots', [tmpz(1) tmpz(1)+100 tmpz(end)], ...
             'leftvalue', 0, ...
             'rightvalue', 75, ...
             'leftslope', 0, ...
             'increasing','on', ...
             'plot','off' ...
             );
% Alternative:
% iodz_lo = findin(odz_lo,zgrid)-5;
%            'rightvalue', [], ... 
%            'knots', [tmpz(1) tmpz(end)], ...

 slm = slmengine(tmpz,tmpv,slm_presc);

 comp.o2_int(iodz_lo:end) = slmeval(tmpz,slm);

 if (1)
    figure
    plot(comp.o2,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.o2_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2) 
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NO2: Split into 3 regions:
 % From inspection of O2 and NO2:
 odz_up = 80;
 odz_lo = 380;
 iodz_up = findin(odz_up,zgrid)-0;
 iodz_lo = findin(odz_lo,zgrid);

 comp.no2_int = nan(size(comp.no2)); 

% (1) oxic upper:
 tmpz = zgrid;
 tmpv = comp.no2;
 
 no2_up = nanmean(comp.no2(1:iodz_up));
 comp.no2_int(1:iodz_up) = no2_up;

 % (3) oxic lower:
 tmpz = zgrid(iodz_lo:end);
 tmpv = comp.no2(iodz_lo:end);

 slm_presc = slmset('degree',3, ...
             'minvalue',0, ... 
             'knots', 3, ...
             'decreasing','on', ...
             'concaveup', 'on', ...
             'plot','off' ...
             );

 slm = slmengine(tmpz,tmpv,slm_presc);

 comp.no2_int(iodz_lo:end) = slmeval(tmpz,slm);

 % Get constraints for intermediate layer: 
 % concentrations, slopes
 no2_lo = slmeval(zgrid(iodz_lo),slm);
 no2_lo_slo = slmeval(zgrid(iodz_lo),slm,1);

% (2) anoxic layer:
 tmpz = zgrid(iodz_up:iodz_lo);
 tmpv = comp.no2(iodz_up:iodz_lo);

 slm_presc = slmset('degree',3, ...
             'knots', [zgrid(iodz_up) zgrid(iodz_up)+50 225 zgrid(iodz_lo)], ...
             'leftvalue', no2_up, ...
             'rightvalue', no2_lo, ...
             'leftslope', 0, ...
             'rightslope', no2_lo_slo, ...
             'minvalue',0, ... 
             'plot','off' ...
             );
            %'knots', 4, ...

 slm = slmengine(tmpz,tmpv,slm_presc);

 comp.no2_int(iodz_up+1:iodz_lo-1) = slmeval(tmpz(2:end-1),slm);

 if (1)
    figure
    plot(comp.no2,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.no2_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2) 
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % PO4
 comp.po4_int = nan(size(comp.po4));

 tmpz = zgrid(1:end);
 tmpv = comp.po4(1:end);

 slm_presc = slmset('degree',3, ...
             'knots', 11, ...
             'leftvalue', comp.po4(1), ...
             'rightslope', 0, ...
             'plot','off' ...
             );
 slm = slmengine(tmpz,tmpv,slm_presc);

 comp.po4_int = slmeval(tmpz,slm);

 if (1)
    figure
    plot(comp.po4,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.po4_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2)
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NO3
 if (0)
    % Straight inteprolation of NO# following PO4
    comp.no3_int = nan(size(comp.no3));
   
    tmpz = zgrid(1:end);
    tmpv = comp.no3(1:end);
   
    slm_presc = slmset('degree',3, ...
                'knots', 11, ...
                'leftvalue', comp.no3(1), ...
                'rightslope', 0, ...
                'plot','off' ...
                );
    slm = slmengine(tmpz,tmpv,slm_presc);
   
    comp.no3_int = slmeval(tmpz,slm);
    comp.nstar_int = 16*comp.po4_int - comp.no3_int;
 else
    % PREFERRED:
    % Interpolates N* and recosntruct NO3
    comp.nstar_int = nan(size(comp.nstar));
   
    tmpz = zgrid(1:end);
    tmpv = comp.nstar(1:end);
   
    slm_presc = slmset('degree',3, ...
                'knots', 12, ...
                'leftvalue', comp.nstar(1), ...
                'rightslope', 0, ...
                'constantregion', [700 tmpz(end)], ...
                'plot','off' ...
                );
    slm = slmengine(tmpz,tmpv,slm_presc);
   
    comp.nstar_int = slmeval(tmpz,slm);
    comp.no3_int = 16*comp.po4_int - comp.nstar_int;
 end

 if (1)
    figure
    plot(comp.no3,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.no3_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2)
    %%%% Check N*
    figure
    plot(comp.nstar,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.nstar_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2)
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % N2O
 if (1)
 % Straight inteprolation with a bunch of best knots
    comp.n2o_int = nan(size(comp.n2o));
   
    tmpz = zgrid(1:end);
    tmpv = comp.n2o(1:end);
   
    slm_presc = slmset('degree',3, ...
                'knots', [30 85 150 300 410 540 665  tmpz(end)], ...
                'leftvalue', comp.n2o(1), ...
                'rightslope', 0, ...
                'plot','off' ...
                );
    slm = slmengine(tmpz,tmpv,slm_presc);
   
    comp.n2o_int = slmeval(tmpz,slm);
 else
 %split domains into upper peak and lower
    odz_up = 180;
    iodz_up = findin(odz_up,zgrid);

    % Uppe part - upper peak
    tmpz = zgrid(1:iodz_up);
    tmpv = comp.n2o(1:iodz_up);

    slm_presc = slmset('degree',3, ...
                'knots', 3, ...
                'leftvalue', comp.n2o(1), ...
                'plot','off' ...
                );
               %'knots', [tmpz(1) 120 tmpz(end)], ...

    slm = slmengine(tmpz,tmpv,slm_presc);
    comp.n2o_int(1:iodz_up) = slmeval(tmpz,slm);

    % Get constraints for lower layer: 
    % concentrations, slopes
    n2o_lo = slmeval(zgrid(iodz_up),slm);
    n2o_lo_slo = slmeval(zgrid(iodz_up),slm,1);

    % Lower part below upper peak
    tmpz = zgrid(iodz_up:end);
    tmpv = comp.n2o(iodz_up:end);
 
    slm_presc = slmset('degree',3, ...
                'knots', [tmpz(1) 280 400 540 665  tmpz(end)], ...
                'leftvalue', n2o_lo, ...
                'leftslope', n2o_lo_slo, ...
                'rightslope', 0, ...
                'plot','off' ...
                );
               %'knots', [tmpz(1) 120 tmpz(end)], ...

    slm = slmengine(tmpz,tmpv,slm_presc);
    comp.n2o_int(iodz_up:end) = slmeval(tmpz,slm);

 end

 if (1)
    figure
    plot(comp.n2o,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.n2o_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2)
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % NH4
 % Cleans up NH4 outlier
 comp.nh4(122) = nan;
 comp.nh4_int = nan(size(comp.nh4));

 tmpz = zgrid(1:end);
 tmpv = comp.nh4(1:end);

 slm_presc = slmset('degree',3, ...
             'knots', 14, ...
             'leftvalue', (comp.nh4(2)+comp.nh4(3))/2, ...
             'linearregion', [200 zgrid(end)], ...
             'plot','off' ...
             );

 slm = slmengine(tmpz,tmpv,slm_presc);

 comp.nh4_int = slmeval(tmpz,slm);

 if (1)
    figure
    plot(comp.nh4,-comp.zgrid,'.b-','markersize',10)
    hold on
    plot(comp.nh4_int,-comp.zgrid,'.r-','markersize',10,'linewidth',2)
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if (1)
    % For consistency with S. Yang code, renames inteprolated
    % variables and keeps both data and interp
    comp.o2_data = comp.o2;
    comp.po4_data = comp.po4;
    comp.no3_data = comp.no3;
    comp.no2_data = comp.no2;
    comp.nh4_data = comp.nh4;
    comp.n2o_data = comp.n2o;
    % Swaps interpolated
    comp.o2 = comp.o2_int;
    comp.po4 = comp.po4_int;
    comp.no3 = comp.no3_int;
    comp.no2 = comp.no2_int;
    comp.nh4 = comp.nh4_int;
    comp.n2o = comp.n2o_int;

    compilation_ETSP_gridded = comp;
    save('compilation_ETSP_gridded_Feb232018_interpol.mat','compilation_ETSP_gridded');
 end
