
function run_epgInversionFSE
  close all; clear;

  J = 50;  % The number of echoes
  t1 = 10;
  t2 = 10;
  T1 = 1000;
  T2 = 100;
  N = 100;
  nomAlpha = 2*pi/3;

  alpha0 = pi/3;
  phi0 = 0;
  Ralpha0 = epgMakeR( alpha0, phi0 );
  Q0m = sparse(3,N);  Q0m(3,1)=1;
	Q0 = sparse( Ralpha0 * Q0m );
  %Q0 = epgFseSim(t1,t2,T1,T2,nomAlpha,3,N);

  figure; hold on;
  for i=1:2
    if i==1
      alpha = nomAlpha * ones(J,1);
    else
      alpha = epgInversionFSE( t1, t2, T1, T2, J, nomAlpha, Q0 );
    end

    % Simulate EPG signal
    allQs = zeros(J,3,N);
    Q = Q0;
    for j=1:J
      Q = epgRelax( Q, t1, T1, T2 );
      Q = epgGrad( Q, 1 );
      Ralpha = epgMakeR( alpha(j), pi/2 );
      Q = Ralpha * Q;
      Q = epgRelax( Q, t2, T1, T2 );
      Q = epgGrad( Q, 1 );
      allQs(j,:,:) = Q;
    end

    sig = squeeze( allQs(:,1,1) );
    plotnice( abs(sig) );
  end
  legend( 'Nominal', 'Adjusted' );
end


function Q = epgFseSim(t1,t2,T1,T2,alpha,J,N)
  alpha0 = pi/2;
  phi0 = 0;
  Q = zeros(3,N);
  Q(3,1) = 1;
  Ralpha0 = epgMakeR( alpha0, phi0 );
  Q = Ralpha0 * Q;
  for j=1:J
    Q = epgRelax( Q, t1, T1, T2 );
    Q = epgGrad( Q, 1 );
    Ralpha = epgMakeR( alpha, pi/2 );
    Q = Ralpha * Q;
    Q = epgRelax( Q, t2, T1, T2 );
    Q = epgGrad( Q, 1 );
  end
end


