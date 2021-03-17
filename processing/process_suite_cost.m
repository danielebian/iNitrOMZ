
 nruns = Suite.nruns;

 [tmp1 tmp2] = bgc1d_fc2minimize_evaluate_cost(Suite.Out{1});
 ncvar = length(tmp2);

 % Initializes outputs
 clear cost;
 cost.all = nan(nruns,1);
 cost.var = nan(nruns,ncvar);

 for indv=1:Suite.nparam
    cost.(Suite.params{indv}) = nan(nruns,1);
 end

 cost.poc100 = nan(nruns,1);
 cost.poc0 = nan(nruns,1);
 
 for indr=1:nruns
    [tmp1 tmp2] = bgc1d_fc2minimize_evaluate_cost(Suite.Out{indr});
    cost.all(indr) = tmp1;
    cost.var(indr,:) = tmp2;
    runindex = cell(Suite.nparam,1);
    [runindex{:}] = ind2sub(Suite.dims,indr);
    for ipar = 1:Suite.nparam
       cost.(Suite.params{ipar})(indr) = Suite.AllParam{ipar}(runindex{ipar});
    end
    cost.poc100(indr) =  Suite.Out{indr}.poc_flux(findin(-100,Suite.Out{1}.zgrid));
    cost.poc0(indr) =  Suite.Out{indr}.poc_flux(1);
 end
 
 % If needed unfolds cost to have same size of structure
 vars = fieldnames(cost);
 nvar = length(vars);
 nsize = Suite.dims;
 for indv=1:nvar
   tmp = cost.(vars{indv});
   osize = size(tmp);
   cost2.(vars{indv}) = reshape(tmp,[nsize osize(2)]);
 end

 % To get indeces of best cost
 [cmin ibest] = min(cost.all);
 bestindex = cell(Suite.nparam,1);
 [bestindex{:}] = ind2sub(Suite.dims,ibest);

 if (1)
 % list first "nmax" best runs
 nmax = 25;
 [scost sind] = sort(cost.all,'ascend');
    for indi=1:nmax
       bestindex = cell(Suite.nparam,1);
       [bestindex{:}] = ind2sub(Suite.dims,sind(indi)); 
       disp([num2str(indi) ' (' num2str(sind(indi)) ') : ' num2str(scost(indi)) ' - ' num2str([bestindex{:}])]); 
    end
 end
 
