% Driver script for solving the 1D advection equations

Globals1D;

% Order of polymomials used for approximation
N = 8;

% Generate simple mesh
K = 20;          %Number of subdomains
xleft = 0.0;
xright = 2*pi;
[Nv, VX, K, EToV] = MeshGen1D(xleft,xright,K);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = sin(x);
%plotting inital condition
    xv = reshape(x,1,[]);
    uv = reshape(u,1,[]);
    plot(xv,uv);
    xlim([xleft xright]);
    ylim([-1.1 1.1]);
    drawnow;

% Solve Problem
FinalTime = 2;
%[u] = Advec1D(u,FinalTime);

%time execution
time = 0;
% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; 
dt = CFL/(2*pi)*xmin; 
dt = .5*dt;
Nsteps = ceil(FinalTime/dt); 
dt = FinalTime/Nsteps;

% advection speed
%a = 2*pi;
a = 4;
%allocating matrix for errors
sx = size(x);
E = zeros(sx(1)*sx(2),Nsteps);

% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        
        %[rhsu] = AdvecRHS1D(u, timelocal, a); original
        % form field differences at faces
        alpha=1.0;
        du = zeros(Nfp*Nfaces,K);
        du(:) = (u(vmapM)-u(vmapP)).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;

        % impose boundary condition at x=0
        uin = -sin(a*time);
        du (mapI) = (u(vmapI)- uin ).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;
        du (mapO) = 0;

        % compute right hand sides of the semi-discrete PDE
        rhsu = -a*rx.*(Dr*u) + LIFT*(Fscale.*(du));
        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end
    
    %analytic
    uana = sin(x-a*time);
    uanav = reshape(uana,1,[]);
        
    % Increment time
    time = time+dt;
    
    %computing error
    xv = reshape(x,1,[]);
    uv = reshape(u,1,[]);
    error = abs(uv - uanav);
    E(:,tstep) = error';
    
    %plotting
    plot(xv,uv);
    hold on;
    plot(xv,uanav);
    hold off;
    xlim([xleft xright]);
    ylim([-1.1 1.1]);
    drawnow;
    
end












