function [sol sadv sdiff ssms srest] = bgc1d_advection_diff(bgc)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Advection-diffusion module
%
% 	Advects, diffuses tracers and applies 
% 	BGC sources minus sinks and restorings fluxes
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Simon Yang, UCLA, April 2019
% Daniele Bianchi, UCLA, June 2019
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

 % Initialize solutions
 sol = zeros(bgc.nt/bgc.hist,bgc.nvar,bgc.nz);
 sadv = zeros(bgc.nt/bgc.hist,bgc.nvar,bgc.nz);
 sdiff = zeros(bgc.nt/bgc.hist,bgc.nvar,bgc.nz);
 ssms = zeros(bgc.nt/bgc.hist,bgc.nvar,bgc.nz);
 srest = zeros(bgc.nt/bgc.hist,bgc.nvar,bgc.nz);

 poc = zeros(2,bgc.nz);
 o2 = zeros(2,bgc.nz);
 no3 = zeros(2,bgc.nz);
 no2 = zeros(2,bgc.nz);
 nh4 = zeros(2,bgc.nz);
 n2o = zeros(2,bgc.nz);
 n2 = zeros(2,bgc.nz);
 po4 = zeros(2,bgc.nz);

 fpoc_out = zeros(2,bgc.nz+1);

 if bgc.RunIsotopes
    i15no3 = zeros(2,bgc.nz);
    i15no2 = zeros(2,bgc.nz);
    i15nh4 = zeros(2,bgc.nz);
    i15n2oA = zeros(2,bgc.nz);
    i15n2oB = zeros(2,bgc.nz);
 end

 if bgc.FromRestart == 1  
    %initial conditions
    o2(1,:) = bgc.rst(1,:);
    no3(1,:) = bgc.rst(2,:);
    poc(1,:) = bgc.rst(3,:);
    po4(1,:) = bgc.rst(4,:);
    n2o(1,:) = bgc.rst(5,:);
    nh4(1,:) = bgc.rst(6,:);
    no2(1,:) = bgc.rst(7,:);
    n2(1,:) = bgc.rst(8,:);  
    if bgc.RunIsotopes
       i15no3(1,:) = bgc.rst(9,:);
       i15no2(1,:) = bgc.rst(10,:);
       i15nh4(1,:) = bgc.rst(11,:);	
       i15n2oA(1,:) = bgc.rst(12,:);
       i15n2oB(1,:) = bgc.rst(13,:);
    end
 else
    poc(1,:) = linspace(bgc.poc_flux_top/bgc.wsink(1),0.01,bgc.npt+1);
    o2(1,:) = linspace(bgc.o2_top,bgc.o2_bot,bgc.npt+1);
    no3(1,:) = linspace(bgc.no3_top,bgc.no3_bot,bgc.npt+1);
    no2(1,:) = 10^-23;linspace(bgc.no2_top,bgc.no2_bot,bgc.npt+1);
    nh4(1,:) = 10^-23;linspace(bgc.nh4_top,bgc.nh4_bot,bgc.npt+1);
    n2o(1,:) = linspace(bgc.n2o_top,bgc.n2o_bot,bgc.npt+1);
    n2(1,:) = linspace(bgc.n2_top,bgc.n2_bot,bgc.npt+1);
    po4(1,:) = linspace(bgc.po4_top,bgc.po4_bot,bgc.npt+1);
    if bgc.RunIsotopes
       i15no3(1,:) = linspace(bgc.i15no3_top,bgc.i15no3_bot,bgc.npt+1);
       i15no2(1,:) = 10^-23*0.0037;linspace(bgc.i15no2_top,bgc.i15no2_bot,bgc.npt+1);
       i15nh4(1,:) = 10^-23*0.0037;linspace(bgc.i15nh4_top,bgc.i15nh4_bot,bgc.npt+1);
       i15n2oA(1,:) = linspace(bgc.i15n2oA_top,bgc.i15n2oA_bot,bgc.npt+1);
       i15n2oB(1,:) = linspace(bgc.i15n2oB_top,bgc.i15n2oB_bot,bgc.npt+1);
    end
 end

 % dump tracers in a structure "tr" - one by one (avoids eval)
 tr.o2 = o2(1,:);
 tr.no3 = no3(1,:);
 tr.poc = poc(1,:);
 tr.po4 = po4(1,:);
 tr.n2o = n2o(1,:);
 tr.nh4 = nh4(1,:);
 tr.no2 = no2(1,:);
 tr.n2 = n2(1,:);
 if bgc.RunIsotopes
    tr.i15no3 = i15no3(1,:);
    tr.i15no2 = i15no2(1,:);
    tr.i15nh4 = i15nh4(1,:);
    tr.i15n2oA = i15n2oA(1,:);
    tr.i15n2oB = i15n2oB(1,:);
 end

 % Get initial SMS
 % Calculate SMS
 % update 15N/N ratios
 sms =  bgc1d_sourcesink(bgc,tr);

 % % % Initialize particulate flux at the top
 fpoc_out(1,1) = bgc.poc_flux_top;

 % % %  Update steady state POC sinking flux
 for z=1:bgc.nz
    % Explicit sinking
    %fpoc_out(1,z+1) = fpoc_out(1,z)*(1.0 - bgc.dz*sms.kpoc(z)/bgc.wsink(z));
    % Implicit sinking
    fpoc_out(1,z+1) = fpoc_out(1,z)/(1.0 + bgc.dz*sms.kpoc(z)/bgc.wsink(z));
    % Updates POC, by mass conservation: remin = divergence of flux
    poc(1,z)=(fpoc_out(1,z)-fpoc_out(1,z+1))/(bgc.dz*sms.kpoc(z));
 end

 % % % % % % % % % % % % % % % % % % 
 % % % % Start time-stepping  % % % %
 % % % % % % % % % % % % % % % % % %
 
 % For advection velocity and diffusion coefficient fixed in time, calculate here
 % terms for the numerical advection-diffusion solver. For time-dependent w and Kv
 % move these terms inside the time loop
 alpha = bgc.wup(2:end-1) * bgc.dt / (2*bgc.dz);
 beta = bgc.dt / (2.*bgc.dz) * (bgc.wup(3:end)-bgc.wup(1:end-2));
 gamma = bgc.Kv(2:end-1) * bgc.dt / (bgc.dz)^2;
 coeff1 = 1+beta-2*gamma;
 coeff2 = alpha+gamma;
 coeff3 = -alpha+gamma;

 for t=1:bgc.nt
    % % % % Now calculate Explicit tracer concentrations
    %%%% Top boundary conditions
    o2(2,1) = bgc.o2_top;
    no3(2,1) = bgc.no3_top;
    no2(2,1) = bgc.no2_top;
    nh4(2,1) = bgc.nh4_top;
    n2o(2,1) = bgc.n2o_top;
    n2(2,1) = bgc.n2_top;
    po4(2,1) = bgc.po4_top;
    %%%% Bottom boundary conditions
    o2(2,end) = bgc.o2_bot;
    no3(2,end) = bgc.no3_bot;
    no2(2,end) = bgc.no2_bot;
    nh4(2,end) = bgc.nh4_bot;
    n2o(2,end) = bgc.n2o_bot;
    n2(2,end) = bgc.n2_bot;
    po4(2,end) = bgc.po4_bot;
    if bgc.RunIsotopes
       %%%% Top boundary conditions
       i15no3(2,1) = bgc.i15no3_top;
       i15no2(2,1) = bgc.i15no2_top;
       i15nh4(2,1) = bgc.i15nh4_top;
       i15n2oA(2,1) = bgc.i15n2oA_top;
       i15n2oB(2,1) = bgc.i15n2oB_top;
       %%%% Bottom boundary conditions
       i15no3(2,end) = bgc.i15no3_bot;
       i15no2(2,end) = bgc.i15no2_bot;
       i15nh4(2,end) = bgc.i15nh4_bot;
       i15n2oA(2,end) = bgc.i15n2oA_bot;
       i15n2oB(2,end) = bgc.i15n2oB_bot;
       % Check for 0 concentrations
       idx=(tr.no3+tr.i15no3==0);
       no3(1,idx)=0;i15no3(1,idx)=0;
       idx=(tr.no2+tr.i15no2==0);
       no2(1,idx)=0;i15no2(1,idx)=0;
       idx=(tr.nh4+tr.i15nh4==0);
       nh4(1,idx)=0;i15nh4(1,idx)=0;
       idx=(tr.n2o+tr.i15n2oA+tr.i15n2oB==0);
       n2o(1,idx)=0;i15n2oA(1,idx)=0;i15n2oB(1,idx)=0;
    end

    %%%% advection and diffusion

    % The code below more compactly and efficiently solve the following equation for all tracers:
    %o2(2,2:end-1) = o2(1,2:end-1) -bgc.wup(2:end-1).*bgc.dt./(2.*-bgc.dz) .* (o2(1,3:end)-o2(1,1:end-2)) - ...
    %                o2(1,2:end-1) .* bgc.dt./(2.*-bgc.dz) .* (bgc.wup(3:end)-bgc.wup(1:end-2)) +  ...
    %                bgc.Kv(2:end-1) .* bgc.dt./(bgc.dz)^2 .* (o2(1,3:end) - 2 .* o2(1,2:end-1) + o2(1,1:end-2)); 

    % Explicitly goes through tracers
    o2(2,2:end-1) = o2(1,2:end-1) .* coeff1 + o2(1,3:end) .* coeff2 + o2(1,1:end-2) .* coeff3;
    no3(2,2:end-1) = no3(1,2:end-1) .* coeff1 + no3(1,3:end) .* coeff2 + no3(1,1:end-2) .* coeff3;
    po4(2,2:end-1) = po4(1,2:end-1) .* coeff1 + po4(1,3:end) .* coeff2 + po4(1,1:end-2) .* coeff3;
    n2o(2,2:end-1) = n2o(1,2:end-1) .* coeff1 + n2o(1,3:end) .* coeff2 + n2o(1,1:end-2) .* coeff3;
    nh4(2,2:end-1) = nh4(1,2:end-1) .* coeff1 + nh4(1,3:end) .* coeff2 + nh4(1,1:end-2) .* coeff3;
    no2(2,2:end-1) = no2(1,2:end-1) .* coeff1 + no2(1,3:end) .* coeff2 + no2(1,1:end-2) .* coeff3;
    n2(2,2:end-1) = n2(1,2:end-1) .* coeff1 + n2(1,3:end) .* coeff2 + n2(1,1:end-2) .* coeff3;
    if bgc.RunIsotopes
       i15no3(2,2:end-1) = i15no3(1,2:end-1) .* coeff1 + i15no3(1,3:end) .* coeff2 + i15no3(1,1:end-2) .* coeff3;
       i15no2(2,2:end-1) = i15no2(1,2:end-1) .* coeff1 + i15no2(1,3:end) .* coeff2 + i15no2(1,1:end-2) .* coeff3;
       i15nh4(2,2:end-1) = i15nh4(1,2:end-1) .* coeff1 + i15nh4(1,3:end) .* coeff2 + i15nh4(1,1:end-2) .* coeff3;
       i15n2oA(2,2:end-1) = i15n2oA(1,2:end-1) .* coeff1 + i15n2oA(1,3:end) .* coeff2 + i15n2oA(1,1:end-2) .* coeff3;
       i15n2oB(2,2:end-1) = i15n2oB(1,2:end-1) .* coeff1 + i15n2oB(1,3:end) .* coeff2 + i15n2oB(1,1:end-2) .* coeff3;
    end

    %%%% Get sources minus Sinks	
    % dump tracers in a structure "tr" - one by one (avoids eval)
    tr.o2 = o2(1,:);
    tr.no3 = no3(1,:);
    tr.poc = poc(1,:);
    tr.po4 = po4(1,:);
    tr.n2o = n2o(1,:);
    tr.nh4 = nh4(1,:);
    tr.no2 = no2(1,:);
    tr.n2 = n2(1,:);
    if bgc.RunIsotopes
       tr.i15no3 = i15no3(1,:);
       tr.i15no2 = i15no2(1,:);
       tr.i15nh4 = i15nh4(1,:);
       tr.i15n2oA = i15n2oA(1,:);
       tr.i15n2oB = i15n2oB(1,:);
    end

    % Calculate SMS
    sms =  bgc1d_sourcesink(bgc,tr);
    
    % % %  Update steady state POC sinking flux
    fpoc_out(2,1) = bgc.poc_flux_top;
    for z=1:bgc.nz
       % Explicit sinking
       %fpoc_out(2,z+1) = fpoc_out(2,z)*(1.0 - bgc.dz*sms.kpoc(z)/bgc.wsink(z));
       % Implicit sinking
       fpoc_out(2,z+1) = fpoc_out(2,z)/(1.0 + bgc.dz*sms.kpoc(z)/bgc.wsink(z));
       % Updates POC, by mass conservation: remin = divergence of flux
       poc(2,z)=(fpoc_out(2,z)-fpoc_out(2,z+1))/(bgc.dz*sms.kpoc(z));
    end       

    %%%% Do sources minus sinks
    o2(2,2:end-1) = o2(2,2:end-1) + sms.o2(2:end-1).*bgc.dt;
    no3(2,2:end-1) = no3(2,2:end-1) + sms.no3(2:end-1).*bgc.dt;
    no2(2,2:end-1) = no2(2,2:end-1) + sms.no2(2:end-1).*bgc.dt;
    nh4(2,2:end-1) = nh4(2,2:end-1) + sms.nh4(2:end-1).*bgc.dt;
    n2o(2,2:end-1) = n2o(2,2:end-1) + sms.n2o(2:end-1).*bgc.dt;
    n2(2,2:end-1) = n2(2,2:end-1) + sms.n2(2:end-1).*bgc.dt;
    po4(2,2:end-1) = po4(2,2:end-1) + sms.po4(2:end-1).*bgc.dt;
    if bgc.RunIsotopes
       i15no3(2,2:end-1) = i15no3(2,2:end-1) + sms.i15no3(2:end-1).*bgc.dt;
       i15no2(2,2:end-1) = i15no2(2,2:end-1) + sms.i15no2(2:end-1).*bgc.dt;
       i15nh4(2,2:end-1) = i15nh4(2,2:end-1) + sms.i15nh4(2:end-1).*bgc.dt;
       i15n2oA(2,2:end-1) = i15n2oA(2,2:end-1) + sms.i15n2oA(2:end-1).*bgc.dt;
       i15n2oB(2,2:end-1) = i15n2oB(2,2:end-1) + sms.i15n2oB(2:end-1).*bgc.dt;
    end
  
    %%%% Do restoring    
    restoring = bgc1d_restoring(bgc,tr);
    o2(2,2:end-1) = o2(2,2:end-1) + restoring.o2(2:end-1).*bgc.dt;
    no3(2,2:end-1) = no3(2,2:end-1) + restoring.no3(2:end-1).*bgc.dt;
    no2(2,2:end-1) = no2(2,2:end-1) + restoring.no2(2:end-1).*bgc.dt;
    nh4(2,2:end-1) = nh4(2,2:end-1) + restoring.nh4(2:end-1).*bgc.dt;
    n2o(2,2:end-1) = n2o(2,2:end-1) + restoring.n2o(2:end-1).*bgc.dt;
    n2(2,2:end-1) = n2(2,2:end-1) + restoring.n2(2:end-1).*bgc.dt;
    po4(2,2:end-1) = po4(2,2:end-1) + restoring.po4(2:end-1).*bgc.dt;    
    if bgc.RunIsotopes
       i15no3(2,2:end-1) = i15no3(2,2:end-1) + restoring.i15no3(2:end-1).*bgc.dt;
       i15no2(2,2:end-1) = i15no2(2,2:end-1) + restoring.i15no2(2:end-1).*bgc.dt;
       i15nh4(2,2:end-1) = i15nh4(2,2:end-1) + restoring.i15nh4(2:end-1).*bgc.dt;
       i15n2oA(2,2:end-1) = i15n2oA(2,2:end-1) + restoring.i15n2oA(2:end-1).*bgc.dt;
       i15n2oB(2,2:end-1) = i15n2oB(2,2:end-1) + restoring.i15n2oB(2:end-1).*bgc.dt;
    end
      
    %%%% old equals new  
    o2(1,:) = o2(2,:);
    no3(1,:) = no3(2,:);
    no2(1,:) = no2(2,:);
    nh4(1,:) = nh4(2,:);
    n2o(1,:) = n2o(2,:);
    n2(1,:) = n2(2,:);
    po4(1,:) = po4(2,:);
    poc(1,:) = poc(2,:);  
    if bgc.RunIsotopes
       i15no3(1,:) = i15no3(2,:);
       i15no2(1,:) = i15no2(2,:);
       i15nh4(1,:) = i15nh4(2,:);
       i15n2oA(1,:) = i15n2oA(2,:);
       i15n2oB(1,:) = i15n2oB(2,:);
    end
  
    %%%% Save history files and diagnostics  
    if mod(t,bgc.hist) == 0
       if bgc.hist_verbose
          disp(['Saving step #' num2str(t/bgc.hist) '/' num2str(bgc.nhist)]);
       end
       % Saving tracer field
       sol((t/bgc.hist),1,:) = o2(1,:);
       sol((t/bgc.hist),2,:) = no3(1,:);
       sol((t/bgc.hist),3,:) = poc(1,:);
       sol((t/bgc.hist),4,:) = po4(1,:);
       sol((t/bgc.hist),5,:) = n2o(1,:);
       sol((t/bgc.hist),6,:) = nh4(1,:);
       sol((t/bgc.hist),7,:) = no2(1,:);
       sol((t/bgc.hist),8,:) = n2(1,:);
       if bgc.RunIsotopes
          sol((t/bgc.hist),9,:) = i15no3(1,:);
          sol((t/bgc.hist),10,:) = i15no2(1,:);
          sol((t/bgc.hist),11,:) = i15nh4(1,:);
          sol((t/bgc.hist),12,:) = i15n2oA(1,:);
          sol((t/bgc.hist),13,:) = i15n2oB(1,:);
       end
    
       %Save fluxes (bgc.flux_diag == 1)
       if bgc.flux_diag == 1

          % Save advection terms
          sadv((t/bgc.hist),1,2:end-1) = alpha .* (o2(1,3:end) - o2(1,1:end-2)) + beta .* o2(1,2:end-1);
          sadv((t/bgc.hist),2,2:end-1) = alpha .* (no3(1,3:end) - no3(1,1:end-2)) + beta .* no3(1,2:end-1);
          sadv((t/bgc.hist),3,2:end-1) = alpha .* (po4(1,3:end) - po4(1,1:end-2)) + beta .* po4(1,2:end-1);
          sadv((t/bgc.hist),4,2:end-1) = alpha .* (n2o(1,3:end) - n2o(1,1:end-2)) + beta .* n2o(1,2:end-1);
          sadv((t/bgc.hist),6,2:end-1) = alpha .* (nh4(1,3:end) - nh4(1,1:end-2)) + beta .* nh4(1,2:end-1);
          sadv((t/bgc.hist),7,2:end-1) = alpha .* (no2(1,3:end) - no2(1,1:end-2)) + beta .* no2(1,2:end-1);
          sadv((t/bgc.hist),8,2:end-1) = alpha .* (n2(1,3:end) - n2(1,1:end-2)) + beta .* n2(1,2:end-1);
          if bgc.RunIsotopes
             sadv((t/bgc.hist),9,2:end-1) = alpha .* (i15no3(1,3:end) - i15no3(1,1:end-2)) + beta .* i15no3(1,2:end-1);
             sadv((t/bgc.hist),10,2:end-1) = alpha .* (i15no2(1,3:end) - i15no2(1,1:end-2)) + beta .* i15no2(1,2:end-1);
             sadv((t/bgc.hist),11,2:end-1) = alpha .* (i15nh4(1,3:end) - i15nh4(1,1:end-2)) + beta .* i15nh4(1,2:end-1);
             sadv((t/bgc.hist),12,2:end-1) = alpha .* (i15n2oA(1,3:end) - i15n2oA(1,1:end-2)) + beta .* i15n2oA(1,2:end-1);
             sadv((t/bgc.hist),13,2:end-1) = alpha .* (i15n2oB(1,3:end) - i15n2oB(1,1:end-2)) + beta .* i15n2oB(1,2:end-1);
          end
          % Save diffusion terms
          sdiff((t/bgc.hist),1,2:end-1) = gamma .* (o2(1,3:end) - 2 * o2(1,2:end-1) + o2(1,1:end-2));
          sdiff((t/bgc.hist),2,2:end-1) = gamma .* (no3(1,3:end) - 2 * no3(1,2:end-1) + no3(1,1:end-2));
          sdiff((t/bgc.hist),3,2:end-1) = gamma .* (po4(1,3:end) - 2 * po4(1,2:end-1) + po4(1,1:end-2));
          sdiff((t/bgc.hist),4,2:end-1) = gamma .* (n2o(1,3:end) - 2 * n2o(1,2:end-1) + n2o(1,1:end-2));
          sdiff((t/bgc.hist),6,2:end-1) = gamma .* (nh4(1,3:end) - 2 * nh4(1,2:end-1) + nh4(1,1:end-2));
          sdiff((t/bgc.hist),7,2:end-1) = gamma .* (no2(1,3:end) - 2 * no2(1,2:end-1) + no2(1,1:end-2));
          sdiff((t/bgc.hist),8,2:end-1) = gamma .* (n2(1,3:end) - 2 * n2(1,2:end-1) + n2(1,1:end-2));
          if bgc.RunIsotopes
             sdiff((t/bgc.hist),9,2:end-1) = gamma .* (i15no3(1,3:end) - 2 * i15no3(1,2:end-1) + i15no3(1,1:end-2));
             sdiff((t/bgc.hist),10,2:end-1) = gamma .* (i15nh4(1,3:end) - 2 * i15nh4(1,2:end-1) + i15nh4(1,1:end-2));
             sdiff((t/bgc.hist),11,2:end-1) = gamma .* (i15nh4(1,3:end) - 2 * i15nh4(1,2:end-1) + i15nh4(1,1:end-2));
             sdiff((t/bgc.hist),12,2:end-1) = gamma .* (i15n2oA(1,3:end) - 2 * i15n2oA(1,2:end-1) + i15n2oA(1,1:end-2));
             sdiff((t/bgc.hist),13,2:end-1) = gamma .* (i15n2oB(1,3:end) - 2 * i15n2oB(1,2:end-1) + i15n2oB(1,1:end-2));
          end

          % Save SMS term
          ssms((t/bgc.hist),1,2:end-1) = sms.o2(2:end-1);
          ssms((t/bgc.hist),2,2:end-1) = sms.no3(2:end-1);
          ssms((t/bgc.hist),4,2:end-1) = sms.po4(2:end-1);
          ssms((t/bgc.hist),5,2:end-1) = sms.n2o(2:end-1);
          ssms((t/bgc.hist),6,2:end-1) = sms.nh4(2:end-1);
          ssms((t/bgc.hist),7,2:end-1) = sms.no2(2:end-1);
          ssms((t/bgc.hist),8,2:end-1) = sms.n2(2:end-1);
          if bgc.RunIsotopes
             ssms((t/bgc.hist),9,2:end-1) = sms.i15no3(2:end-1);
             ssms((t/bgc.hist),10,2:end-1) = sms.i15no2(2:end-1);
             ssms((t/bgc.hist),11,2:end-1) = sms.i15nh4(2:end-1);
             ssms((t/bgc.hist),12,2:end-1) = sms.i15n2oA(2:end-1);
             ssms((t/bgc.hist),13,2:end-1) = sms.i15n2oB(2:end-1);
          end
        
          % Save restoring term
          srest((t/bgc.hist),1,2:end-1) = restoring.o2(2:end-1);
          srest((t/bgc.hist),2,2:end-1) = restoring.no3(2:end-1);
          srest((t/bgc.hist),4,2:end-1) = restoring.po4(2:end-1);
          srest((t/bgc.hist),5,2:end-1) = restoring.n2o(2:end-1);
          srest((t/bgc.hist),6,2:end-1) = restoring.nh4(2:end-1);
          srest((t/bgc.hist),7,2:end-1) = restoring.no2(2:end-1);
          srest((t/bgc.hist),8,2:end-1) = restoring.n2(2:end-1);
          if bgc.RunIsotopes
             srest((t/bgc.hist),9,2:end-1) = restoring.i15no3(2:end-1);
             srest((t/bgc.hist),10,2:end-1) = restoring.i15no2(2:end-1);
             srest((t/bgc.hist),11,2:end-1) = restoring.i15nh4(2:end-1);
             srest((t/bgc.hist),12,2:end-1) = restoring.i15n2oA(2:end-1);
             srest((t/bgc.hist),13,2:end-1) = restoring.i15n2oB(2:end-1);
          end  % bgc.RunIsotopes
       end  % bgc.flux_diag
    end  % mod(t,bgc.hist)
 end  % t

 % Save restart (bgc.SaveRestart == 1)
 if bgc.SaveRestart == 1
    rst = squeeze(sol(end,:,:));
    endtime = num2str(bgc.nt*bgc.dt/3600/24/365,'%5.1f');
    save([bgc.root, '/restart/', bgc.RunName,'_restart_',endtime,'.mat'],'rst');
 end	

