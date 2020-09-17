function mcost = bgc1d_fc2minimize_cmaes_parallel(ParStart,FunArg)

 % Unfolds arguments
 ParNames = FunArg.ParNames;
 ParMin = FunArg.ParMin;
 ParNorm = FunArg.ParNorm;

 % Re-builds non-normalized parameters:
%ParVal = ParMin + ParStart .* ParNorm;
 ParVal = bsxfun(@plus,ParMin,bsxfun(@times,ParStart,ParNorm));

 % Allows for parallelization of cost function, accepting ParStart input of size NpxNm
 % Np: # parameters of the problem
 % Nm: # simultaneous calls

 [Np Nm] = size(ParStart);

 mcost = nan(1,Nm);

 parfor indm=1:Nm

    % Initialize model:
    bgc = bgc1d_initialize;
    bgc.hist_verbose = false;
   
    % Change parameters with those selected in Scenario_child(ichr).chromosome
    for indp=1:length(ParNames)
       bgc = bgc1d_change_input(bgc,ParNames{indp},ParVal(indp,indm));
    end

    % % % % % % % % % % % %
    % DERIVED PARAMETERS  %
    % % % % % % % % % % % %
    % Here imposes the following constraint:
    % Constraints: KDen1 + KDen2 + KDen3 = remin
%    remin = 0.08/86400;
%    slowFactor =3.3; % BEC
%    bgc.KDen1 = max(0, remin/slowFactor - bgc.KDen2 - bgc.KDen3);
    % % % % % % % % % % % %
   
    % Update dependent parameters
    bgc = bgc1d_initialize_DepParam(bgc);
    if bgc.RunIsotopes
       bgc = bgc1d_initIso_Dep_params(bgc);
    end
   
    % % % % % % %
    % Run model %
    % % % % % % %
   %disp(['Running bgc1d_iteration']);
    [bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
    bgc.sol = squeeze(bgc.sol_time(end,:,:));
    
    % Get the cost by calling the cost function
    % Note on the options:
    % (0) will run without output
    % (1) DON'T USE (will create a figure on screen)
    % (2) Will save out a jpeg evaluation of the cost for each run
    tcost = bgc1d_fc2minimize_evaluate_cost(bgc,0);
   %tcost = bgc1d_fc2minimize_evaluate_cost(bgc,2);
    
    % Fill in parallel array of costs
    mcost(1,indm) = tcost;
    disp(['bgc1d_iteration - cost : ' num2str(tcost)]);

 end

