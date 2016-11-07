
function run_epgInversionFSE
  close all; clear;

  t1 = 100;
  t2 = 115;

  T1 = 1000;
  T2 = 800;

  alpha = epgInversionFSE( t1, t2, T1, T2 );


end

