function tcost = bgc1d_fc2minimize_cmaes(ParStart,FunArg)

 % Unfolds arguments
 ParNames = FunArg.ParNames;
 ParMin = FunArg.ParMin;
 ParNorm = FunArg.ParNorm;

 % Original list of parameters tested by Simon:
 % {'KO2Rem','KNO2No' 'KO2Den1','KO2Den2' 'KO2Den3','KNO3Den1','KDen1','KDen2', ...
 %  'KDen3','KNO2Den2','KN2ODen3','KAx','KNH4Ax','KNO2Ax','KO2Ax'};

 % Re-builds non-normalized parameters:
 ParVal = ParMin + ParStart .* ParNorm;

 % Initialize model:
 bgc = bgc1d_initialize;
 bgc.hist_verbose = false;

 % Change parameters with those selected in Scenario_child(ichr).chromosome
 for indp=1:length(ParNames)
    bgc = bgc1d_change_input(bgc,ParNames{indp},ParVal(indp));
 end

    % % % % % % % % % % % %
    % DERIVED PARAMETERS  %
    % % % % % % % % % % % %
    % Here imposes the following constraint:
    % Constraints: KDen1 + KDen2 + KDen3 = remin
    remin = 0.08/86400;
    bgc.KDen1 = max(0, remin - bgc.KDen2 - bgc.KDen3);
    % % % % % % % % % % % %

 % Update dependent parameters
 bgc = bgc1d_initialize_DepParam(bgc);
 if bgc.RunIsotopes
    bgc = bgc1d_initIso_Dep_params(bgc);
 end

 % % % % % % %
 % Run model %
 % % % % % % %
 [bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
 bgc.sol = squeeze(bgc.sol_time(end,:,:));

 % Get the cost by calling the cost function
 tcost = bgc1d_fc2minimize_evaluate_cost(bgc);
 disp(['bgc1d_iteration - cost : ' num2str(tcost)]);

