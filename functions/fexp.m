  function fexp = fexp(lvar,klvar);

  % If needed turn on following line to prevent less than zero variable values
  % (note that this is not necessary since the same operation is done in sms)
  %lvar = max(0,lvar);
   fexp = exp(-lvar/klvar);

  end
