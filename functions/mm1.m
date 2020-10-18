  function fmen1 = mmen1(var,Kvar);

  % If needed turn on following line to prevent less than zero variable values
  % (note that this is not necessary since the same operation is done in sms)
  %var = max(0,var);
   fmen1 = var./(var+Kvar);

  end
