
function alpha = epgInversionFSE( t1, t2, T1, T2 )

  J = 1;  % number of echoes in FSE sequence
  F0D = 0.1;   % Desired signal strength at all echoes
  theta = 30 * pi/180;  % max value of alpha
  N = 4;  % largest EGP coefficient to store
  phi = 0;

  alpha0 = findAlpha0( F0D, theta );
  Ralpha0 = epgMakeR( alpha0, phi-pi/2 );  % Meet CPMG conditions
  Q0m = sparse(3,N);  Q0m(3,1)=1;
	Q0 = sparse( Ralpha0 * Q0m );

  E11 = exp( -t1 / T1 );
  E22 = exp( -t2 / T2 );
  E1122 = E11 * E22;
  E1 = exp( -(t1+t2)/T1 );
  E2 = exp( -(t1+t2)/T2 );

  %b = makebVec( t1, t2, T1, T2, J, Q0, F0D );
  %H = makeH( t1, t2, T1, T2, J, Q0, phi );

  megaH = zeros(2*J,J);
  megab = zeros(2*J,1);
  for thisJ=1:J
    b = makebVec( E11, E2, thisJ, Q0, F0D );
    megab(1+(thisJ-1)*2:2+(thisJ-1)*2) = b;
    H = makeH( E11, E22, E1122, E1, E2, thisJ, Q0, phi );
    megaH(1+(thisJ-1)*2:2+(thisJ-1)*2) = H;
  end

  alpha = megaH \ megab;
end


function alpha0 = findAlpha0( F0d, theta, Qin )
  if nargin < 3, Qin=[0;0;1;]; end;

  dAlpha = 1/360*pi;
  bestAlpha0 = 0;
  cost = 99999;
  for alpha = -theta:dAlpha:theta
    Ralpha0 = epgMakeR( alpha, -pi/2 );
    thisQ = Ralpha0 * Qin;
    if abs(thisQ(1,1)-F0d) < cost
      bestAlpha0 = alpha;
      cost = abs(thisQ(1,1)-F0d);
    end
  end
  alpha0 = bestAlpha0;
end


function b = makebVec( E11, E2, M, Q0, F0D )
  QE = sparse( size(Q0,1), size(Q0,2) ); QE(3,1) = 1;

  A = makeA( E11, E2, size(Q0,2) );

  Q0_vec = [ real(Q0(:)); imag(Q0(:)); ];
  QE_vec = [ real(QE(:)); imag(QE(:)); ];
  tmp = A^M * Q0_vec;
  for i=2:M
    tmp = tmp + A^(M-i+1) * QE_vec;
  end

  b = [ F0D; 0;] - [ tmp(1,1); tmp(numel(Q0)+1,1); ];
end


function H = makeH( E11, E22, E1122, E1, E2, M, Q0, phi )
  % t1 is the time until the RF pulse
  % t2 is the time from the RF pulse to the end of the repetition
  % T1 and T2 are the decay constants
  % M is the number of echoes
  % Q0 is the initial state
  % phi is the phase of the RF pulse

  A = makeA( E11, E2, size(Q0,2) );
  B = makeB( E1122, E1, size(Q0,2) );

  QE = sparse( size(Q0,1), size(Q0,2) );  QE(3,1) = 1;
  QR = sparse( size(Q0,1), size(Q0,2) );  QR(1,2) = E22*(1-E11)*(-1i*exp(1i*phi));

  preH = cell(M,M);
  Q0_vec = [ real(Q0(:)); imag(Q0(:)); ];  
  for i=1:M
    preH{1,i} = A^(M-i) * B * A^(i-1) * Q0_vec;
  end

  QE_vec = [ real(QE(:)); imag(QE(:)); ];
  for j=1:M
    for i=j:M
      preH{j,i} = A^(M-i) * B * A^(i-j) * QE_vec;
    end
  end

  QR_vec = [ real(QR(:)); imag(QR(:)); ];
  for i=2:M
    preH{i,i} = preH{i,i} + A^(M-i+1) * QR_vec;
  end

  H = cell(1,M);
  for i=1:M
    H{i} = sparse( size(Q0_vec,1), size(Q0_vec,2) );
    for j=1:i
      H{i} = H{i} + preH{j,i};
    end
  end

  H = cell2mat(H);
  H = [ H(1:M,:); H(numel(Q0)+1,:); ];
end


function A = makeA( E11, E2, N )
  % N is the index of the largest EPG coefficient vector retained
  AFn_real = [ circshift( E2*eye(N), 2, 2 ), zeros(N) ];
  AFn_real(end-1:end,:)=0;
  AFn_imag = [ zeros(N), circshift( E2*eye(N), 2, 2 ) ];
  AFn_imag(end-1:end,:) = 0;

  AFp_real = [ circshift( E2*eye(N), -2, 2 ), zeros(N) ];
  AFp_real(1,:)=0;  AFn_real(1,3)= E2;  % AFn here is not a mistake
  AFp_real(2,:)=0;  AFn_real(2,2)= E2;
  AFp_imag = [ zeros(N), circshift( E2*eye(N), -2, 2 ) ];
  AFp_imag(1,:)=0;  AFp_imag(1,N+3)= -E2;
  AFp_imag(2,:)=0;  AFp_imag(2,N+2)= -E2;

  AZ = [ E11*eye(N), zeros(N); zeros(N), E11*eye(N) ];

  A = sparse( 2*3*N, 2*3*N );
  A(1:3:end,1:3:end) = [ AFp_real; AFp_imag; ];
  A(2:3:end,2:3:end) = [ AFn_real; AFn_imag; ];
  A(3:3:end,3:3:end) = AZ;
end


function B = makeB( E1122, E1, N )
  % N is the index of the largest EPG coefficient vector retained
  BFn_real = [ circshift( E1122*eye(N), 1, 2 ), zeros(N) ];
  BFn_real(end,:)=0;
  BFn_imag = [ zeros(N), circshift( E1122*eye(N), 1, 2 ) ];
  BFn_imag(end,:)=0;

  BFp_real = [ circshift( -E1122*eye(N), -1, 2 ), zeros(N) ];
  BFp_real(1,:)=0;  BFp_real(1,2) = -E1122;
  BFp_imag = [ zeros(N), circshift( -E1122*eye(N), -1, 2 ) ];
  BFp_imag(1,:)=0;  BFp_imag(1,N+2)= E1122;

  BZm_real = [ circshift( E1*eye(N), 1, 2 ), zeros(N) ];
  BZm_real(end,:)=0;
  BZm_imag = [ zeros(N), circshift( E1*eye(N), 1, 2 ) ];
  BZm_imag(end,:)=0;

  BZp_real = [ circshift( -E1*eye(N), -1, 2 ), zeros(N) ];
  BZp_real(1,:)=0;  BZm_real(1,2) = -E1;  % BZm here is not a mistake
  BZp_imag = [ zeros(N), circshift( -E1*eye(N), -1, 2 ) ];
  BZp_imag(1,:)=0;  BZm_imag(2,N+2) = E1;  % BZm here is not a mistake

  B = sparse( 2*3*N, 2*3*N );
  BFp = [ BFp_real; BFp_imag; ];  B(1:3:end,3:3:end) = BFp;
  BFn = [ BFn_real; BFn_imag; ];  B(2:3:end,3:3:end) = BFn;
  BZ1 = [ BZm_real; BZm_imag; ];  B(3:3:end,2:3:end) = BZ1;
  BZ2 = [ BZp_real; BZp_imag; ];  B(3:3:end,1:3:end) = BZ2;
end


function AQ = applyA( Q, E11, E2 )
  AQ = E11 * Q;  % only last row is correct here
  AQ(2,1:end-2) = E2 * Q(2,3:end);
  AQ(2,end-2:end) = 0;
  AQ(1,3:end) = E2 * Q(1,1:end-2);
  AQ(1,2) = E2 * conj( Q(2,2) );
  AQ(1,1) = E2 * conj( Q(1,1) );
end

function BQ = applyB( Q, E1122, E1 )
  BQ = sparse( size(Q,1), size(Q,2) );
  BQ(3,2:end-1) = E1 * ( Q(2,3:end) - Q(1,1:end-2) );
  BQ(3,1) = E1 * ( Q(2,2) - conj( Q(1,2) ) );
  BQ(2,1:end-1) = E1122 * Q(3,2:end);
  BQ(1,2:end) = -E1122 * Q(3,1:end-1);
  BQ(1,1) = -E1122 * conj( Q(3,2) );
end

