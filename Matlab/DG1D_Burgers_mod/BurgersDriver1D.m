% clear;
clc;
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
% Professor Jan S. Hesthaven
%
% PDE:
%           du/dt + du^2/dx = 0  for -1 < x < +1
%         IC:
%           epsilon = 0.1
%           u(x,0) = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0

% clear
% clc
% Globals1D;

filtereo = [2,1,0,2,1,0,2,1,0];
%2 -> Traditional filter
%1 -> Energy based filter
%0 -> No filter
epsilones = [0.0005,0.0005,0.0005,0.00025,0.00025,0.00025,0.000001,0.000001,0.000001];

elementos = [6,9,12,15];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%CONTROL PANEL
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For analyticial computation of uo=-sin(pi*x)
[etai,wi] = JacobiGQ(0,0,2048);
% normuv = zeros(13000,18);
% normuv = normuv*NaN;
N=22;
solut = {zeros(N+1,elementos(1)),zeros(N+1,elementos(1)),zeros(N+1,elementos(2)),...
    zeros(N+1,elementos(2)),zeros(N+1,elementos(3)),zeros(N+1,elementos(3)),...
    zeros(N+1,elementos(4)),zeros(N+1,elementos(4))};


for iii=1:4

% Order of polymomials used for approximation
% N = 22;

%Final time
FinalTime = 0.6;

%dealiasing??
deal = 0.0;     %1.0 for dealiasing, 0.0 for non-dealiasing

%improved integration??
iint = 0.0;     %1.0 for immproved integration, 0.0 for not
%hay un error con iint=1.0, dan mal los flujos viscosos

%init Filter matrix
%filter = ?
%1.0 for filtering, 0.0 for non-filtering
% if filtereo(iii) == 2
%     filter = 1;
% else
    filter = 0;
% end
% if filtereo(iii) == 1
    filtere = 1;
% else
%     filtere = 0;
% end

% t=0 for Non-Unitary, 1 for unitary, 2 for Lobatto Basis.
% t = 0;          
% if t==2
%     Nc=4;           %when using Lobatto Bassis you should put Nc=3 for preserve 2 first basis functions (because keep the B.C. of subdomains)
%     s=1;           % exponential filter degree
% else
%     Nc = 1;
%     s=6;           % exponential filter degree
% end

%viscosity constant
% epsilon = 0.0;      %original
% epsilon = 0.001;      %No need of filtering
% epsilon = 0.0005;      %viscosity constant (need of filtering and works)
% caso de éxito es epsilon=0.0005, Elements=12, N=16, zeros and
% filterdiag(i), con RHSm de test
% epsilon = 0.00025;      %viscosity constant (need of filtering and does not work)
epsilon = 0.000001;
% epsilon = epsilones(iii);

% Generate simple mesh
xL = -1;
xR = 1;
% Elements = 10;
Elements = elementos(iii);
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,Elements);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%END CONTROL PANEL
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Initialize solver and construct grid and metric
StartUp1D;

% Filter matrix F
if filter==1.0
    %last parameter (t) means type,
    Ncf = 1;
    s = N/2;
    [F,L] = Filter1D(N,Ncf,s,0,V,W,r);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set initial conditions
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;    %original Intitial C.
%u = (x<-0.5)*2;     %step initial condition
u = -sin(2*pi*x/(xR-xL));
% u = x;
% uab = 2;
% uaa = 1;
% xao = -0.5;
% u(x<=xao) = uab;
% u(x>xao) = uaa;
uo = u;
ua = u;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Compute time step size
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  xmin = min(abs(x(1,:)-x(2,:)));
  if epsilon <= 0.0001
  CFL = 0.1;
  else
  CFL = 1.0;
  end
  umax = max(max(abs(u)));
  dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));     %original
%   dt = 0.1*CFL*(xmin/umax);
%   dt = CFL*(xmin/umax);
%   dt = 0.0005;
  Nsteps = ceil(FinalTime/dt);
  dt = FinalTime/Nsteps; 

  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Initial Energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Init_energy
  modes = (0:N)';

%%%%%%%%
%limits for integration within analytical
if epsilon > 0.001
    ai = -3;
else
    ai = -2;
end
%%%%%%%%

%
%  Outer time step loop 
%
%   figure ( 1 );
%   shapestr = { '-o','-x' };

for tstep=1:Nsteps

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Going into RK4
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     if time >= 0.4
%         pause
%     end
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        %[rhsu] = BurgersRHS1D(u,epsilon,xL,xR,timelocal);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BurgersRHS1D subroutine
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Define field differences at faces of the K elements.
        %
          du = zeros(Nfp*Nfaces,K); 
          du(:) = u(vmapM) - u(vmapP);
        %
        %  Impose boundary condition at x=0
        %
          %uin =-tanh((xL+0.5-time) / (2*epsilon)) + 1.0;   %original
          uin = 0;                     %Dirichlet = 0
%           uin = uab;                    %Dirichlet for step function
%           uin = u((N+1)*Elements);      %Periodic
          du(mapI) = 2.0*(u(vmapI)-uin);
          %uout=-tanh((xR+0.5-time) / (2*epsilon)) + 1.0;   %original
          uout = 0;                    %Dirichlet = 0
%           uout = uaa;                   %Dirichlet for step function
%           uout = u((N+1)*Elements);     %Periodic
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
%          tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
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
          %dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q );       %original
          
          if deal==1.0
              %umd = invV*u;             %u modal
              
              %u modal and zero-padding
              umd = [invV*u ; zeros(length(rd)-length(r),Elements)];
              
              %going to nodal with more nodes
              ud = Vd*umd;
              
              %computing non-linear term with more nodes
              dfdrd = Drd * ( ud .^ 2 / 2);
              %dfdrd = ud .^ 2;
              
              %settling the non-linear result in modal
              umd = invVd*dfdrd;
              
              %substraction the zero padding part
              umd = umd(1:N+1,:);
              
              %going back to original nodal
              ud = V*umd;

              %computing the discretized term
              dfdr = ud - Dr * (sqrt ( epsilon ) * q );
              %dfdr = Dr * ( (ud ./ 2) - sqrt ( epsilon ) * q );  
          else
              dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q );
          end
          
        %
        %  Compute right hand sides of the semi-discrete PDE
        %
          rhsu = - ( rx .* dfdr - LIFT * ( Fscale .* flux ) );
        
        %solving velocity  
          resu = rk4a(INTRK)*resu + dt*rhsu;
          u = u+rk4b(INTRK)*resu;
        
%           %Ensure continuity at element interfaces
%             for k = 2:K
%                 u(Np,k-1) = (u(1,k) + u(Np,k-1))/2;
%                 u(1,k) = u(Np,k-1);
%                 %uf(Np,k-1) = u(Np,k-1);        %uncomment for not filter the boundaries of the elements
%                 %uf(1,k) = u(1,k);              %uncomment for not filter the boundaries of the elements
%             end
          
        % end BurgersRHS1D subroutine
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % energy calculations
          Calc_energy             
        %end energy calculations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
%     %computing analytical solution
    ua = analitica(x,time+dt,epsilon,etai,wi,ai);      %-sin(pi*x)
%     ua = analitica2(x,time+dt,uab,uaa,xao);          %Inviscid Burgers with step
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % energy storage
% Stop condition (sometimes useful)
%     if time > 0.1
%         pause
%     end
    Store_energy
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Filtering    
    if filter == 1.0
        u = filtering(u,F);
        %u = F*u;
    end
    
%    % storing norm
%    normuv(tstep,iii*3-2) = time;
%    normuv(tstep,iii*3-1) = norm(abs(u-ua));
%    normuv(tstep,iii*3) = norm(abs(u-ua),'Inf');
%     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Advancing time  
    time = time + dt;

%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Plotting
%     %
%     %  Display solution on every 10th time step.
%     %
%     if ( rem ( tstep, 10 ) == 0 )
% 
%       % Plotting analytical and numerical solutions  
% %       subplot(1,2,1)
%       %ploting analytical
%       plot(x,ua,'k')
%       hold on
%       for ii = 1 : Elements
% %         plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
% %           'Markersize', 1, 'LineWidth', 2.0, 'Color','b');
%         plot ( x(:,ii), u(:,ii), shapestr{1+rem(ii,2)}, ...
%           'Markersize', 1, 'LineWidth', 1.5);
%         hold all
%       end
%       grid ( 'on' );
%       axis ( [ xL, xR, uBott, uUp ] );
%       
%       drawnow;
%       hold off
%       
% %      plot ( xv, ustv, 'Markersize', 1, 'LineWidth', 1.5, 'Color','r');
%       
% %       % Plotting energy by frequency
% %       subplot(1,2,2)
% %       plot(modes,Etom);
% %       xlim([-0.1 N+0.1]);
% %       ylim([0 0.2]);
% %       legend('1','2','3','4','5','6','7','8','9','10','11','12')
% %       drawnow
%       
%     end
    %pause(0.05)
    
    %ploting
% %     xv = reshape(x,1,[]);
% %     uv = reshape(u,1,[]);
%     subplot(1,2,1)
%     plot(x,u);
%     xlim([xL xR]);
%     ylim([uBott uUp]);
%     drawnow

end

solut{1,iii*2-1} = x;
solut{1,iii*2} = u;

end
%Postprocessing

%end Burgers1D subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


