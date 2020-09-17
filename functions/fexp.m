  function oxygen_function = f2o2(o2,ko2);

   o2pos = max(0,o2);
   oxygen_function = exp(-o2pos/ko2);

  end
