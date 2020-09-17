 function bgc = bgc1d_initbgc_params(bgc)

 d2s = 86400;
 %%%%%%%%% Stochiometry %%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Organic matter form: C_aH_bO_cN_dP
 % Stoichiometric ratios: C:H:O:N:P = a:b:c:d:1
 % Anderson and Sarmiento 1994 stochiometry
 bgc.stoch_a = 106.0;
 bgc.stoch_b = 175.0;
 bgc.stoch_c = 42.0;
 bgc.stoch_d = 16.0;

 %%%%%%% Ammonification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.Krem = 0.08/d2s ;              % 0.08    % Max. remineralization rate (1/s)
 bgc.KO2Rem  = 0.046571922641370;   % 4       % Half sat. constant for respiration  (mmolO2/m3) - Martens-Habbena 2009

 %%%%%% Ammonium oxidationn %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Ammox: NH4 --> NO2
 bgc.KAo = 0.123585996046118/d2s;       % 0.045    Max. Ammonium oxidation rate (1/s) - Bristow 2017
 bgc.KNH4Ao  = 0.128718397274300;       % 0.1       % Half sat. constant for nh4 (mmolN/m3) - Peng 2016
 bgc.KO2Ao = 0.333;                     % 0.333+-0.130  % Half sat. constant for Ammonium oxidation (mmolO2/m3) - Bristow 2017

 %%%%%%% Nitrite oxidationn %%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Nitrox: NO2 --> NO3
 bgc.KNo = 0.05/d2s;                  % 0.256    Max. Nitrite oxidation rate (1/s) - Bristow 2017
 bgc.KNO2No = 1.0;                    % Don't know (mmolN/m3)
 bgc.KO2No = 0.778;                   % Half sat. constant of NO2 for Nitrite oxidation (mmolO2/m3) - Bristow 2017

 %%%%%%% Denitrification %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Denitrif1: NO3 --> NO2
 bgc.KDen1 = 0.024105944781783/d2s;     % Max. denitrif1 rate (1/s)
 bgc.KO2Den1 = 5.0;                     % O2 poisoning constant for denitrif1 (mmolO2/m3)
 bgc.KNO3Den1 = 0.422400394303500;      % Half sat. constant of NO3 for denitrif1 (mmolNO3/m3)

 % Denitrif2: NO2 --> N2O
 bgc.KDen2 = 0.008/d2s;                 % Max. denitrif2 rate (1/s)
 bgc.KO2Den2 = 2.0;                     % O2 poisoning constant for denitrif2 (mmolO2/m3)
 bgc.KNO2Den2 = 0.01;                   % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3)

 % Denitrif3: N2O --> N2
 bgc.KDen3 = 0.055693889684983/d2s;     % Max. denitrif3 rate (1/s)
 bgc.KO2Den3 = 0.611650781313400;       % O2 poisoning constant for denitrif3 (mmolO2/m3)
 bgc.KN2ODen3 = 0.2;                    % Half sat. constant of N2O for denitrification3 (mmolNO3/m3)

 %%%%%%%%%% Anammox %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.KAx = 0.395422078044298/d2s;       % Max. Anaerobic Ammonium oxidation rate (1/s) - Bristow 2017
 bgc.KNH4Ax  = 0.230901186513400;       % Half sat. constant of NH4 for anammox (mmolNH4/m3)
 bgc.KNO2Ax  = 0.1;                     % Half sat. constant of NO2 for anammox (mmolNO2/m3)
 bgc.KO2Ax = 0.01;                      % 1.0     %

 %%%%% N2O prod via ammox %%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Parameters for calculation of N2O yields during Ammox and
 % nitrifier-denitrification (see n2o_yield.m).
 
 % Choose paramterization:
 % 'Ji': Ji et al 2018
%bgc.n2o_yield = 'Ji';
 bgc.n2o_yield = 'Yang';

 switch bgc.n2o_yield
 case 'Ji'
    % Ji et al 2018
    bgc.Ji_a = 0.2;
    bgc.Ji_b  = 0.08;
 case 'Yang'
    % S. Yang March 15 2019 optimized results 
    bgc.Ji_a = 0.45;
    bgc.Ji_b  = 0.19;
 otherwise
    error(['N2O yield ' bgc.n2o_yield  ' case not found']);
 end


