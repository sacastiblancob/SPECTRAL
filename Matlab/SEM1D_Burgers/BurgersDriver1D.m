% Driver script for solving the 1D burgers equations
%
% PDE:
%           du/dt + du^2/dx = 0  for -1 < x < +1
%         IC:
%           epsilon = 0.1
%           u(x,0) = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0


Globals1D;

% Order of polymomials used for approximation
N = 8;

% Generate simple mesh
xL = -1.0;
xR = 1.0;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
epsilon = 0.1;
u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;

% Solve Problem
FinalTime = 1.5;
[u] = Burgers1D ( u, epsilon, xL, xR, FinalTime );



