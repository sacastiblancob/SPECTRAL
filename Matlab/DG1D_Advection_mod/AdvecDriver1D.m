clear; clc;
% Driver script for solving the 1D advection equations

Globals1D;

% Order of polymomials used for approximation
N = 16;
Elements = 6;

% Generate simple mesh
xL = -1.0;
xR = 1.0;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,Elements);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
% u = sin(x);
% sin function
% u = sin(2*pi*x)+1;

% Gaussian bump
% u = exp(-((x-(-0.0)).^2)/(2*(0.1^2)));

%Step function
u = zeros(size(x));
% u(ceil((N+1)/2+2*N):ceil((N+1)/2+3*N)) = 1.0;
u(10:20) = 1.0;

% Solve Problem
FinalTime = 2;
%[u] = Advec1D(u,FinalTime);

time = 0.0;

%
%  Runge-Kutta residual storage.
%
resu = zeros ( Np, K );
%
%  Compute time step size.
%
xmin = min ( abs ( x(1,:) - x(2,:) ) );
CFL = 0.75;
dt = CFL / ( 2 * pi ) * xmin;
dt = 0.5 * dt;
Nsteps = ceil ( FinalTime / dt );
dt = FinalTime / Nsteps;
%
%  Advection speed.
%

%a = 2.0 * pi;
a = 2.0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Initial Energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Init_energy

%
%  Outer time step loop.
%
figure ( 1 );
shapestr = { '-o','-x' };

for tstep = 1 : Nsteps
%     if tstep==2
%       pause
%     end
    for INTRK = 1 : 5
          timelocal = time + rk4c ( INTRK ) * dt;
          %[rhsu] = AdvecRHS1D ( u, timelocal, a );

          alpha = 1.0;
          du = zeros(Nfp*Nfaces,K);
          %du(:) = (u(vmapM)-u(vmapP)).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;     %original
          %du(:) = u(vmapM).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2 - u(vmapP).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;
          du(:) = u(vmapM).*(a*nx(:))/2 - u(vmapM).*((1-alpha)*abs(a*nx(:)))/2 - u(vmapP).*(a*nx(:))/2 + ...
              u(vmapP).*((1-alpha)*abs(a*nx(:)))/2;
        %
        %  Impose boundary condition at x=0
        %
          %ORIGINAL
           uin = u((N+1)*Elements);    %periodic boundary condition for a > 0
%           uin = -sin(2*pi*a*time);   %analitical solution for impose boundary condition
%           uin = 0;
          du(mapI) = (u(vmapI)- uin ).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;   %original
%           du(mapI) = u(vmapI).*(a*nx(mapI))/2 - u(vmapI)*((1-alpha)*abs(a*nx(mapI)))/2 - uin.*(a*nx(mapI))/2 + ...
%              uin.*((1-alpha)*abs(a*nx(mapI)))/2;

          du(mapO) = 0; %original

        %
        %  Compute right hand sides of the semi-discrete PDE
        %
          rhsu = -a*rx.*(Dr*u) + LIFT*(Fscale.*(du));
          %rhsu = ((1-x.^2).^5 + 1).*(Dr*u) + LIFT*(Fscale.*(du));
          %rhsu = u.*(Dr*u) + LIFT*(Fscale.*(du));  %%jajaja

          resu = rk4a ( INTRK ) * resu + dt * rhsu;
          u = u + rk4b ( INTRK ) * resu;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % energy calculations
          Calc_energy             
        %end energy calculations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % energy storage        
    Store_energy
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %
    %  Increment time.
    %
    time = time + dt;
    %
    %  Display solution on every 10th time step.
    %
    if ( rem ( tstep, 10 ) == 0 )
      for i = 1 : Elements
        plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
          'Markersize', 8, 'LineWidth', 2 );
        hold all
      end
      grid ( 'on' );
      axis ( [ xL, xR, -3.0, 3.0 ] );
      drawnow;
      hold off
    end

end

Econ = EEt + nfEEt;
%Econ = EEt + dEEt + nfEEt;

Postprocessing

%end Advection1D subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





