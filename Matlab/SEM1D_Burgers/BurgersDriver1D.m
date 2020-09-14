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

% epsilon = sqrt(viscosity)
%epsilon = 0.01;
%nu = epsilon^2;
nu = 0.05;
epsilon = sqrt(nu);

%init Filter matrix
%filter = ?
filter = 0.0;   %1.0 for filtering, 0.0 for non-filtering
t = 2;          % 0 for Non-Unitary, 1 for unitary, 2 for Lobatto Basis.
Nc=3;           %when using Lobatto Bassis you should put Nc=3 for preserve 2 first basis functions (because keep the B.C. of subdomains)
if t==2
    s=15;           % exponential filter degree
else
    s=10;           % exponential filter degree
end

% Generate simple mesh
xL = -1.0;
xR = 1.0;
elements = 10;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,elements);

% Initialize solver and construct grid and metric
StartUp1D;

if filter==1.0
    %last parameter (t) means type,
    F = Filter1D(N,Nc,s,t,V,W,r);
end

% Set initial conditions
%u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;
%uv = exp(-xv.^2/sqrt(0.01));
uv = -sin(pi*xv);

% Solve Problem
FinalTime = 1.5;
%[u] = Burgers1D ( u, epsilon, xL, xR, FinalTime );

time = 0.0;
%
%  Runge-Kutta residual storage  
%
  %resu = zeros(Np,K);
  resu = zeros(Ns,1);
%
%  Compute time step size
%
  xmin = min(abs(x(1,:)-x(2,:)));
  CFL = 0.25;
  umax = max(max(abs(uv)));
  dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));
  Nsteps = ceil(FinalTime/dt);
  dt = FinalTime/Nsteps; 
%
%  Outer time step loop 
%
  vel = 1.1;
  for tstep=1:Nsteps

    for INTRK = 1:5
      timelocal = time + rk4c(INTRK)*dt;
      %[rhsu] = BurgersRHS1D(u,epsilon,xL,xR,timelocal);
      %uv(1) = exp(-(xv(1)-vel*timelocal).^2/sqrt(0.01));
      uv(1) = 0;
      uv(Ns) = 0;
      %uv(Ns) = exp(-(xv(Ns)-vel*timelocal).^2/sqrt(0.01));
      rhsu = -1.0.*AG*(uv.^2);
      if nu ~= 0
        rhsu = (nu*SG)\rhsu;
      end
      resu = rk4a(INTRK)*resu + dt*rhsu;
      uv = uv+rk4b(INTRK)*resu;
    end
    time = time + dt;
    
    if filter==1.0
        uv = filtering(uv,F,K,Np);
    end
        
    
    
    %ploting
%     xv = reshape(x,1,[]);
%     uv = reshape(u,1,[]);
    plot(xv,uv);
    xlim([xL xR]);
    ylim([-1.1 1.1]);
    drawnow
    
  end


