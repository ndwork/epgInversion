
function alphaOut = epgInversionFSE( t1, t2, T1, T2, J, alpha, Q0 )
  % J is the number of echoes in FSE sequence

  F0mD = 1.0;   % Desired signal strength at all echoes
  N = 100;  % largest EGP coefficient to store
  %theta = 30 * pi/180;  % max value of alpha
  theta = pi/2;

  E11 = exp( -t1 / T1 );
  E12 = exp( -t1 / T2 );
  E21 = exp( -t2 / T1 );
  E22 = exp( -t2 / T2 );
  
  Ra = epgMakeR( alpha, pi/2 );  % phi0=0, phi=pi/2

  A = makeA( E11, E12, E21, E22, size(Q0,2), Ra );  
  B = makeB( E11, E12, E21, E22, size(Q0,2), Ra );

  Qtilde = sparse( size(Q0,1), size(Q0,2) );
  Qtilde(1,2) = E22*(1-E11)*sin(alpha);

  megaH = zeros(2*J,J);
  megab = zeros(2*J,1);
  for thisJ=1:J
    disp([ 'Working on echo ', num2str(thisJ), ' of ', num2str(J) ]);
    b = makebVec( A, B, thisJ, Q0, Qtilde, alpha, F0mD );
    megab(1+(thisJ-1)*2:2+(thisJ-1)*2) = b;
    H = makeH( A, B, thisJ, Q0, Qtilde );
    megaH(1+(thisJ-1)*2:2+(thisJ-1)*2,1:size(H,2)) = [ H(1,:); H(numel(Q0)+1,:); ];
  end
  megab = real( megab );
  megaH = real( megaH );

  useLsqLin = 1;
  if useLsqLin == 1,
    lb = -theta * ones( J, 1 );
    ub = theta * ones( J, 1 );
    dAlpha = lsqlin( megaH, megab, [], [], [], [], lb, ub );
  else
    opts = tfocs;
    opts.alg = 'N83';
    opts.maxIts = 200;
    opts.printEvery = 1;
    opts.tol = 1d-8;
    x0 = zeros( J, 1 );
    [dAlpha,tfocsOptOut] = tfocs( smooth_quad, { megaH, -megab }, ...
      proj_linf(theta), x0, opts );
    %[dAlpha,tfocsOptOut] = tfocs( smooth_huber(0.1), { megaH, -megab }, ...
    %  proj_linf(theta), x0, opts );
  end

  alphaOut = alpha + dAlpha;
end


function bVec = makebVec( A, B, M, Q0, Qtilde, alpha, F0D )
  Q0_vec = [ real(Q0(:)); imag(Q0(:)); ];
  Qtilde_vec = [ real(Qtilde(:)); imag(Qtilde(:)); ];

  tmp = A^M * Q0_vec;
  for i=1:M
    tmp = tmp - A^(M-i)*B*A^i * Q0_vec * alpha;
  end

  for j=1:M
    tmp = tmp + A^(M-j) * Qtilde_vec;

    for i=j+1:M
      tmp = tmp - A^(M-i) * B * A^i * Qtilde_vec;
    end
  end

  bVec = [ 0; F0D;] - [ tmp(1,1); tmp(numel(Q0)+1,1); ];
end


function H = makeH( A, B, M, Q0, Qtilde )
  % M is the number of echoes
  % Q0 is the initial state

  preH = cell(M,M);
  Q0_vec = [ real(Q0(:)); imag(Q0(:)); ];
  for i=1:M
    preH{1,i} = A^(M-i) * B * A^(i-1) * Q0_vec;
  end

  Qtilde_vec = [ real(Qtilde(:)); imag(Qtilde(:)); ];
  for j=2:M
    for i=j:M
      preH{j,i} = A^(M-i) * B * A^(i-j) * Qtilde_vec;
    end
  end

  H = cell(1,M);
  for i=1:M
    H{i} = sparse( size(Q0_vec,1), size(Q0_vec,2) );
    for j=1:i
      H{i} = H{i} + preH{j,i};
    end
  end

  H = cell2mat(H);
end


function A = makeA( E11, E12, E21, E22, N, Ra )
  % N is the index of the largest EPG coefficient vector retained

  % Notation for matrices: A<out><in>
  AFpFn_real = E22 * Ra(1,2) * E12 * [ eye(N), zeros(N) ];
  AFpFn_imag = E22 * Ra(1,2) * E12 * [ zeros(N), eye(N) ];

  AFpFp_real = E22 * Ra(1,1) * E12 * [ circshift(eye(N),-2,2), zeros(N) ];
  AFpFp_real(1:2,:)=0;
  AFpFp_imag = E22 * Ra(1,1) * E12 * [ zeros(N), circshift(eye(N),-2,2) ];
  AFpFp_imag(1:2,:)=0;

  AFpFn_real(1,3) = AFpFn_real(1,3) + E22 * Ra(1,1) * E12;
  AFpFn_real(2,2) = AFpFn_real(2,2) + E22 * Ra(1,1) * E12;
  AFpFn_imag(1,3) = AFpFn_imag(1,3) - E22 * Ra(1,1) * E12;
  AFpFn_imag(2,2) = AFpFn_imag(2,2) - E22 * Ra(1,1) * E12;

  AFpZ_real = E22 * Ra(1,3) * E11 * [ circshift(eye(N),-1,2), zeros(N) ];
  AFpZ_real(1,:)=0;  AFpZ_real(1,2)=1;
  AFpZ_imag = E22 * Ra(1,3) * E11 * [ zeros(N), circshift(eye(N),-1,2) ];
  AFpZ_imag(1,:)=0;  AFpZ_imag(1,2)= -1;

  AFnFp_real = E22 * Ra(2,1) * E12 * [ eye(N), zeros(N) ];
  AFnFp_imag = E22 * Ra(2,1) * E12 * [ zeros(N), eye(N) ];

  AFnFn_real = E22 * Ra(2,2) * E12 * [ circshift(eye(N),2,2), zeros(N) ];
  AFnFn_real(end-1:end,:)=0;
  AFnFn_imag = E22 * Ra(2,2) * E12 * [ zeros(N), circshift(eye(N),2,2) ];
  AFnFn_imag(end-1:end,:)=0;

  AFnZ_real = E22 * Ra(2,3) * E11 * [ circshift(eye(N),1,2), zeros(N) ];
  AFnZ_real(end,:)=0;
  AFnZ_imag = E22 * Ra(2,3) * E11 * [ zeros(N), circshift(eye(N),1,2) ];
  AFnZ_imag(end,:)=0;

  AZFn_real = E21 * Ra(3,2) * E12 * [ circshift(eye(N),1,2), zeros(N) ];
  AZFn_real(end,:)=0;
  AZFn_imag = E21 * Ra(3,2) * E12 * [ zeros(N), circshift(eye(N),1,2) ];
  AZFn_imag(end,:)=0;

  AZFp_real = E21 * Ra(3,1) * E12 * [ circshift(eye(N),-1,2), zeros(N) ];
  AZFp_real(1,:)=0;
  AZFp_imag = E21 * Ra(3,1) * E12 * [ zeros(N), circshift(eye(N),-1,2) ];
  AZFp_imag(1,:)=0;

  AZFn_real(1,2) = AZFn_real(1,2) + E21 * Ra(3,1) * E12;
  AZFn_imag(1,2) = AZFn_imag(1,2) - E21 * Ra(3,1) * E12;

  AZZ_real = E21 * Ra(3,3) * E11 * [ eye(N), zeros(N) ];
  AZZ_imag = E21 * Ra(3,3) * E11 * [ zeros(N), eye(N) ];

  A = sparse( 2*3*N, 2*3*N );
  A(1:3:end,1:3:end) = A(1:3:end,1:3:end) + [ AFpFp_real; AFpFp_imag; ];
  A(1:3:end,2:3:end) = A(1:3:end,2:3:end) + [ AFpFn_real; AFpFn_imag; ];
  A(1:3:end,3:3:end) = A(1:3:end,3:3:end) + [ AFpZ_real; AFpZ_imag; ];
  A(2:3:end,1:3:end) = A(2:3:end,1:3:end) + [ AFnFp_real; AFnFp_imag; ];
  A(2:3:end,2:3:end) = A(2:3:end,2:3:end) + [ AFnFn_real; AFnFn_imag; ];
  A(2:3:end,3:3:end) = A(2:3:end,3:3:end) + [ AFnZ_real; AFnZ_imag; ];
  A(3:3:end,1:3:end) = A(3:3:end,1:3:end) + [ AZFp_real; AZFp_imag; ];
  A(3:3:end,2:3:end) = A(3:3:end,2:3:end) + [ AZFn_real; AZFn_imag; ];
  A(3:3:end,3:3:end) = A(3:3:end,3:3:end) + [ AZZ_real; AZZ_imag; ];
end


function B = makeB( E11, E12, E21, E22, N, Ra )
  % N is the index of the largest EPG coefficient vector retained

  % Notation for matrices: B<out><in>
  BFpFn_real = E22 * Ra(1,3) * -E12/2 * [ eye(N), zeros(N) ];
  BFpFn_imag = E22 * Ra(1,3) * -E12/2 * [ zeros(N), eye(N) ];

  BFpFp_real = E22 * Ra(1,3) * -E12/2 * [ circshift(eye(N),-2,2), zeros(N) ];
  BFpFp_real(1:2,:)=0;
  BFpFp_imag = E22 * Ra(1,3) * -E12/2 * [ zeros(N), circshift(eye(N),-2,2) ];
  BFpFp_imag(1:2,:)=0;
  BFpFn_real(1,3) = BFpFn_real(1,3) + E22 * Ra(1,3) * -E12/2;
  BFpFn_imag(1,3) = BFpFn_imag(1,3) - E22 * Ra(1,3) * -E12/2;
  BFpFn_real(2,2) = BFpFn_real(2,2) + E22 * Ra(1,3) * -E12/2;
  BFpFn_imag(2,2) = BFpFn_imag(2,2) - E22 * Ra(1,3) * -E12/2;

  BFpZ_real = E22 * E11 * ( Ra(1,1) + Ra(1,2) ) *  [ circshift(eye(N),-1,2), zeros(N) ];
  BFpZ_real(1,:)=0;  BFpZ_real(1,2) = E22 * E11 * ( Ra(1,1) + Ra(1,2) );
  BFpZ_imag = E22 * E11 * ( Ra(1,1) + Ra(1,2) ) * [ zeros(N), circshift(eye(N),-1,2) ];
  BFpZ_imag(1,:)=0;  BFpZ_imag(1,2) = -E22 * E11 * ( Ra(1,1) + Ra(1,2) );

  BFnFp_real = E22 * Ra(2,3) * -E12/2 * [ eye(N), zeros(N) ];
  BFnFp_imag = E22 * Ra(2,3) * -E12/2 * [ zeros(N), eye(N) ];

  BFnFn_real = E22 * Ra(2,3) * -E12 * [ circshift(eye(N),2,2), zeros(N) ];
  BFnFn_real(end-1:end,:) = 0;
  BFnFn_imag = E22 * Ra(2,3) * -E12 * [ zeros(N), circshift(eye(N),2,2) ];
  BFnFn_imag(end-1:end,:) = 0;

  BFnZ_real = E22 * E11 * ( Ra(2,1) + Ra(2,2) ) * [ circshift(eye(N),1,2), zeros(N) ];
  BFnZ_real(end,:) = 0;
  BFnZ_imag = E22 * E11 * ( Ra(2,1) + Ra(2,2) ) * [ zeros(N), circshift(eye(N),1,2) ];
  BFnZ_imag(end,:) = 0;

  BZFn_real = E21 * Ra(3,3) * -E12/2 * [ circshift(eye(N),1,2), zeros(N) ];
  BZFn_real(end,:) = 0;
  BZFn_imag = E21 * Ra(3,3) * -E12/2 * [ zeros(N), circshift(eye(N),1,2) ];
  BZFn_imag(end,:) = 0;

  BZFp_real = E21 * Ra(3,3) * -E12/2 * [ circshift(eye(N),-1,2), eye(N) ];
  BZFp_real(1,:) = 0;
  BZFp_imag = E21 * Ra(3,3) * -E12/2 * [ eye(N), circshift(eye(N),-1,2) ];
  BZFp_imag(1,:) = 0;

  BZFn_real(1,2) = BZFn_real(1,2) + E21 * Ra(3,3) * -E12/2;
  BZFn_imag(1,2) = BZFn_imag(1,2) - E21 * Ra(3,3) * -E12/2;

  BZZ_real = E21 * E11 * ( Ra(3,1) + Ra(3,2) ) * [ eye(N), zeros(N) ];
  BZZ_imag = E21 * E11 * ( Ra(3,1) + Ra(3,2) ) * [ zeros(N), eye(N) ];

  B = sparse( 2*3*N, 2*3*N );
  B(1:3:end,1:3:end) = B(1:3:end,1:3:end) + [ BFpFp_real; BFpFp_imag; ];
  B(1:3:end,2:3:end) = B(1:3:end,2:3:end) + [ BFpFn_real; BFpFn_imag; ];
  B(1:3:end,3:3:end) = B(1:3:end,3:3:end) + [ BFpZ_real; BFpZ_imag; ];
  B(2:3:end,1:3:end) = B(2:3:end,1:3:end) + [ BFnFp_real; BFnFp_imag; ];
  B(2:3:end,2:3:end) = B(2:3:end,2:3:end) + [ BFnFn_real; BFnFn_imag; ];
  B(2:3:end,2:3:end) = B(2:3:end,2:3:end) + [ BFnZ_real; BFnZ_imag; ];
  B(3:3:end,1:3:end) = B(3:3:end,1:3:end) + [ BZFp_real; BZFp_imag; ];
  B(3:3:end,2:3:end) = B(3:3:end,2:3:end) + [ BZFn_real; BZFn_imag; ];
  B(3:3:end,3:3:end) = B(3:3:end,3:3:end) + [ BZZ_real; BZZ_imag; ];
end

