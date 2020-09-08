% Driver script for solving the 1D Maxwell’s equations
Globals1D;

% Polynomial order used for approximation
N = 6;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-2.0,2.0,80);

% Initialize solver and construct grid and metric
StartUp1D;

% Set up material parameters
eps1 = [ones(1,K/2), 2*ones(1,K/2)];
mu1 = ones(1,K);
epsilon = ones(Np,1)*eps1; mu = ones(Np,1)*mu1;

% Set initial conditions
E = sin(pi*x).*(x<0); H = zeros(Np,K);

% Solve Problem
FinalTime = 10;
%[E,H] = Maxwell1D(E,H,epsilon,mu,FinalTime);

time = 0;

% Runge-Kutta residual storage
resE = zeros(Np,K); resH = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=1.0;
dt = CFL*xmin;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        %[rhsE, rhsH] = MaxwellRHS1D(E,H,eps,mu);
                
        % Compute impedance
        Zimp = sqrt(mu./eps);

        % Define field differences at faces
        dE = zeros(Nfp*Nfaces,K);
        dE(:) = E(vmapM)-E(vmapP);
        dH = zeros(Nfp*Nfaces,K);
        dH(:) = H(vmapM)-H(vmapP);
        Zimpm = zeros(Nfp*Nfaces,K);
        Zimpm(:) = Zimp(vmapM);
        Zimpp = zeros(Nfp*Nfaces,K);
        Zimpp(:) = Zimp(vmapP);
        Yimpm = zeros(Nfp*Nfaces,K);
        Yimpm(:) = 1./Zimpm(:);
        Yimpp = zeros(Nfp*Nfaces,K);
        Yimpp(:) = 1./Zimpp(:);

        % Homogeneous boundary conditions, Ez=0
        Ebc = -E(vmapB);
        dE (mapB) = E(vmapB) - Ebc;
        Hbc = H(vmapB);
        dH (mapB) = H(vmapB) - Hbc;

        % evaluate upwind fluxes
        fluxE = 1./(Zimpm + Zimpp).*(nx.*Zimpp.*dH - dE);
        fluxH = 1./(Yimpm + Yimpp).*(nx.*Yimpp.*dE - dH);

        % compute right hand sides of the PDE’s
        rhsE = (-rx.*(Dr*H) + LIFT*(Fscale.*fluxE))./eps;
        rhsH = (-rx.*(Dr*E) + LIFT*(Fscale.*fluxH))./mu;
        
        resE = rk4a(INTRK)*resE + dt*rhsE;
        resH = rk4a(INTRK)*resH + dt*rhsH;
        E = E+rk4b(INTRK)*resE;
        H = H+rk4b(INTRK)*resH;
    end
    % Increment time
    time = time+dt;
end



