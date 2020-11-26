%INITIAL INPUTS alpha, beta and N
%   N = polynomial order
%   if alpha=beta=0.0 --> Legendre Polynomials
%   if alpha=beta=1.0 --> Chebyshev Polynomials
%
% and remind that if you have an polynomial order N you have N+1 nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.0;
beta = 0.0;
p = 6;   %polynomial grade

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GLL nodes, vanderomende and GLL weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nodes
x = JacobiGL(alpha,beta,p);
N = length(x);

%vandermonde matrix
B = Vandermondes(x,alpha,beta);
B3 = legendreM(x,N);
%B is the unitary ones
%B3 is just orthogonal

%GLL weights
w = 2./(N*(N-1)*(B3(:,N).^2));
W = diag(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proofs functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computing u
un = x.^2;
um = B\un;

%prooving equality
%here works
left = B'*W*B*um;
rig = B'*W*un;

%here does not work
left3 = B3'*W*B3*um;
rig3 = B3'*W*un;

%computing the integral between -1 and 1
INT = un'*w;

%unitary
U = B'*W*B;
C = inv(U);
U3 = B3'*W*B3;
C3 = inv(U3);





