% Driver script for solving the 1D burgers equations
%
% PDE:
%           du/dt + 0.5*d(u^2)/dx = 0  for -1 < x < +1
%         IC:
%           epsilon = 0.1
%           u(x,0) = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0


Globals1D;

%%%GENERAL CONFIGURATION
configuration;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,Elements);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
epsilon = 0.005;
u = Initial_condition(x, epsilon);
uo = u;
ua = u;

% Solve Problem
FinalTime = 3/pi;
Burgers1D
% [u] = Burgers1D ( u, epsilon, xL, xR, FinalTime );

if ebf==1
uv_ebf = uv;

if filter_ebf == 1
Error_exp = Error;
Err2_t_exp = Err2_t;
Errinf_t_exp = Errinf_t;
Errord_t_exp = Errord_t;
else
Error_van = Error;
Err2_t_van = Err2_t;
Errinf_t_van = Errinf_t;
Errord_t_van = Errord_t;
end

end

if cf==1
uv_cf = uv;
Err2_t_cf_4 = Err2_t;
Errinf_t_cf_4 = Errinf_t;
Errord_t_cf_4 = Errord_t;
Error_cf_4 = Error;
end

if svv==1
uv_svv = uv;
Err2_t_svv = Err2_t;
Errinf_t_svv = Errinf_t;
Errord_t_svv = Errord_t;
Error_svv = Error;
end




