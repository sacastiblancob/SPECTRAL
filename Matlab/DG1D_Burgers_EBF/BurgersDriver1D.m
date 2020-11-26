%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Licensing:
%
%    Permission to use this software for noncommercial
%    research and educational purposes is hereby granted
%    without fee.  Redistribution, sale, or incorporation
%    of this software into a commercial product is prohibited.
%
%    THE AUTHORS OR PUBLISHER DISCLAIMS ANY AND ALL WARRANTIES
%    WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
%    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
%    PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHORS OR
%    THE PUBLISHER BE LIABLE FOR ANY SPECIAL, INDIRECT OR
%    CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
%    RESULTING FROM LOSS OF USE, DATA OR PROFITS.
%
%  Modified:
%
%    20 September 2018
%
%  Author:
%
%    Original version by Jan Hesthaven, Tim Warburton.
%    Some modifications by John Burkardt.
%
%  Reference:
%
%    Jan Hesthaven, Tim Warburton,
%    Nodal Discontinuous Galerkin Methods: 
%    Algorithms, Analysis, and Applications,
%    Springer, 2007,
%    ISBN: 978-0387720654.
%

% Driver script for solving the 1D burgers equations
%
% PDE:
%           du/dt + du^2/dx = 0  for -1 < x < +1
%         IC:
%           epsilon = 0.1
%           u(x,0) = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Globals1D;

% Order of polymomials used for approximation
N = 8;

% Generate simple mesh
xL = -1.0;
xR = 1.0;
K = 10;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,K);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
epsilon = 0.0;
u = -sin(2*pi*x/(xR-xL));
% u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;

%Down and Up limits
uBott = min(min(u)) - 0.1;
uUp = max(max(u)) + 0.1;

% Solve Problem
FinalTime = 1.5;

  time = 0.0;
%
%  Runge-Kutta residual storage  
%
  resu = zeros(Np,K); 
%
%  Compute time step size
%
  xmin = min(abs(x(1,:)-x(2,:)));
  CFL = 0.25;
  umax = max(max(abs(u)));
  dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));
  Nsteps = ceil(FinalTime/dt);
  dt = FinalTime/Nsteps; 

%
%  Energy initialization
%
  EBF_init  
  
%
%  Outer time step loop 
%
  figure ( 1 );
  shapestr = { '-o','-x' };

  for tstep=1:Nsteps
    %if(tstep==170)
    %    return
    %end
    for INTRK = 1:5
        
        %Local time for RK4
        timelocal = time + rk4c(INTRK)*dt;
      
        %
        %  Define field differences at faces of the K elements.
        %
          du = zeros(Nfp*Nfaces,K); 
          du(:) = u(vmapM) - u(vmapP);
        %
        %  Impose boundary condition at x=0
        %
          % uin =-tanh((xL+0.5-time) / (2*epsilon)) + 1.0; 
          uin = 0.0;
          du(mapI) = 2.0*(u(vmapI)-uin);
          % uout=-tanh((xR+0.5-time) / (2*epsilon)) + 1.0; 
          uout = 0.0;
          du(mapO) = 2.0*(u(vmapO) - uout);
        %
        %  Compute Q.
        %
          q = sqrt ( epsilon ) * ( rx .* ( Dr * u ) - LIFT*(Fscale.*(nx.*du/2.0)));
        %
        %  Compute jumps DQ at each element face.
        %
          dq = zeros(Nfaces,K); 
          dq(:) = ( q(vmapM) - q(vmapP) ) / 2.0;
        %
        %  Impose boundary conditions on the jumps.
        %  Here, Dirichlet conditions are used.
        %
          dq(mapI) = 0.0; 
          dq(mapO) = 0.0;
        %
        %  Evaluate DU2, nonlinear flux at faces of the K elements.
        %
          du2 = zeros(Nfp*Nfaces,K); 
          du2(:) = ( u(vmapM).^2 - u(vmapP).^2 ) / 2.0;
        % 
        %  Impose boundary conditions on the fluxes.
        %
          du2(mapI) = ( u(vmapI).^2 - uin.^2 ); 
          du2(mapO) = ( u(vmapO).^2 - uout.^2 );
        %
        %  Compute flux
        %
          maxvel = max ( max ( abs ( u ) ) );
        %
        %  Penalty scaling -- See Chapter 7.2
        %
        %  tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
        %
          tau = 0.0;
        %
        %  Flux term
        %
          flux = nx .* ( du2 / 2.0 - sqrt ( epsilon ) * dq ) ...
            - maxvel / 2.0 .* du - sqrt ( epsilon ) * tau .* du;
        %
        %  DFDD, local derivatives of field
        %
          dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q ); 
        %
        %  Compute right hand sides of the semi-discrete PDE
        %
          rhsu = - ( rx .* dfdr - LIFT * ( Fscale .* flux ) );

        %Solving U with RK4
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
        
        %
        % Energy Calculations and Filtering
        %
          EBF_calc
        
    end
    time = time + dt;
    
%     %ploting
%     xv = reshape(x,1,[]);
%     uv = reshape(u,1,[]);
%     plot(xv,uv);
%     xlim([xL xR]);
%     ylim([0 2.1]);
%     drawnow
    
    if ( rem ( tstep, 10 ) == 0 )

      % Plotting analytical and numerical solutions  
%      subplot(1,2,1)
%       %ploting analytical
%       plot(x,ua,'k')
%       hold on
      for i = 1 : K
        plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
          'Markersize', 1, 'LineWidth', 2 );
        hold all
      end
      grid ( 'on' );
      axis ( [ xL, xR, uBott, uUp ] );
      
      drawnow;
      hold off
      pause(0.1)
%       % Plotting energy by frequency
%       subplot(1,2,2)
%       plot(modes,Etom);
%       xlim([-0.1 N+0.1]);
%       ylim([0 0.2]);
%       legend('1','2','3','4','5','6','7','8','9','10','11','12')
%       drawnow
      
    end
    
  end


