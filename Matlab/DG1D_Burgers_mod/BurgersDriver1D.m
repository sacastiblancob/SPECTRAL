% Driver script for solving the 1D burgers equations
%
% PDE:
%           du/dt + du^2/dx = 0  for -1 < x < +1
%         IC:
%           epsilon = 0.1
%           u(x,0) = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0


Globals1D;

% Order of polymomials used for approximation
N = 12;

% Generate simple mesh
xL = -1;
xR = 1;
Elements = 12;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,Elements);

% Initialize solver and construct grid and metric
StartUp1D;

%init Filter matrix
%filter = ?
filter = 1.0;

if filter==1.0
    Nc=3;   %when using Lobatto Bassis you should put Nc=3 for preserve 2 first basis functions (because keep the B.C. of subdomains)
    s=36;
    t = 0;  % 0 for Non-Unitary, 1 for unitary, 2 for Lobatto Basis.
    %last parameter (t) means type,
    F = Filter1D(N,Nc,s,t);
end

% Set initial conditions
%epsilon = 0.05;      %original
%epsilon = 0.000001;      %viscosity constant
epsilon = 0.0;
%u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;    %original Intitial C.
%u = (x<-0.5)*2;     %step initial condition
u = -sin(2*pi*x/(xR-xL));
EE = zeros(1,K);
EEdf = zeros(1,K);
EEnf = zeros(1,K);

%Down and Up limits
uBott = min(min(u)) - 0.1;
uUp = max(max(u)) + 0.1;

% Solve Problem
FinalTime = 1.2;
%[u] = Burgers1D ( u, epsilon, xL, xR, FinalTime );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Burgers1D subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%init time
time = 0.0;
%
%  Runge-Kutta residual storage  
%
resu = zeros(Np,K);
rese = zeros(1,K);      %result vector for energy disipation due to viscosity
resedf = zeros(1,K);    %result vector for energy disipation flux over bound.
resenf = zeros(1,K);    %result vector for energy flux due to non-linear term
%
%  Compute time step size
%
xmin = min(abs(x(1,:)-x(2,:)));
CFL = 0.25;
umax = max(max(abs(u)));
%dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));     %original
%dt = 0.1*CFL*(xmin/umax);
dt = CFL*(xmin/umax);
%dt = 3.1776e-04;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps; 

%
%  Initial Energy
%

%Initialization
E = zeros(Nsteps+1,1);      %Energy of whole domain
Ev = zeros(Nsteps+1,1);     %Energy dissipated by viscosity in whole domain
Evf = zeros(Nsteps+1,1);    %Fluxes of energy due to viscosity in whole domain
Enlf = zeros(Nsteps+1,1);   %Non-linear fluxes of energy in whole domain
dEEt = zeros(Nsteps+1,K);   %energy dissipation by element due to viscosity
dfEEt = zeros(Nsteps+1,K);  %energy flux due to dissipation
nfEEt = zeros(Nsteps+1,K);  %energy flux due to non-linear term
EEt = zeros(Nsteps+1,K);    %energy by element
%Econ = zeros(Nsteps+1,K);   %conservation?
T = zeros(Nsteps+1,1);

%Initial Energy
Eto = u.^2;
EEt(1,:) = w*(Eto.*J);
E(1) = sum(EEt(1,:));
T(1) = 0.0;


%
%  Outer time step loop 
%
  figure ( 1 );
  shapestr = { '-o','-x' };

for tstep=1:Nsteps

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
          %uin = u((N+1)*Elements);      %Periodic
          du(mapI) = 2.0*(u(vmapI)-uin);
          %uout=-tanh((xR+0.5-time) / (2*epsilon)) + 1.0;   %original
          uout = 0;                    %Dirichlet = 0
          %uout = u((N+1)*Elements);     %Periodic
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
        
        % end BurgersRHS1D subroutine
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % energy calculations
        
        %disipation over boundaries
        %Fscale = 1/jacobiano en las fronteras de los elementos
        deru = Dr*u./J;
        dflux = 2*epsilon.*u([1 N+1],:).*deru([1 N+1],:);
        dflux = (dflux(2,:) - dflux(1,:));
        
        %non-linear flux over boundaries
        nflux = 2*(u([1 N+1],:).^3)./3;
        nflux = -(nflux(2,:) - nflux(1,:));
        
        %dissipation term
        edis = -(2*epsilon*w*(((deru).^2).*J));
        
      %solving energy
      %RHSE of energy equation integrated over a domain
      %rhse = edis; %dflux + nflux + edis;
      
      %dissipation due to viscous term
      rese = rk4a(INTRK)*rese + dt*edis;
      EE = EE + rk4b(INTRK)*rese;
      
      %viscous fluxes term
      resedf = rk4a(INTRK)*resedf + dt*dflux;
      EEdf = EEdf + rk4b(INTRK)*resedf;
      
      %non-linear fluxes term
      resenf = rk4a(INTRK)*resenf + dt*nflux;
      EEnf = EEnf + rk4b(INTRK)*resenf;      
            
      %end energy calculations
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      %solving velocity  
      resu = rk4a(INTRK)*resu + dt*rhsu;
      u = u+rk4b(INTRK)*resu;
    end
    
%     %computing analytical solution
%     ua = analitica(x,time+dt,epsilon);
        
    % computing energy before filtering
    Et = u.^2;
    %Eele = 0.5*w*(Et.*J);
    EEt(tstep+1,:) = w*(Et.*J);
    E(tstep+1) = sum(EEt(tstep+1,:));
    
%     %Computing dissipate energy due to viscosity
%     %duv = (1./J).*Dr*u;
%     dus = (Dr*u./J).^2;
%     %function of N, Elements, CFL, xL and xR
%     %Evele = 2*0.008059650327027*epsilon*w*(dus.*J);
%     Evele = dt*epsilon*w*(dus.*J);
%     %Evele = epsilon*w*(dus.*J);
    dEEt(tstep+1,:) = -EE;
    dfEEt(tstep+1,:) = -EEdf;
    nfEEt(tstep+1,:) = -EEnf;
    Ev(tstep+1) = -1.0*(sum(EE));
    Evf(tstep+1) = -1.0*(sum(EEdf));
    Enlf(tstep+1) = -1.0*(sum(EEnf));
    
%    Evtot = sum(EE);
%    Ev(tstep+1) = Ev(tstep) + Evtot;
        
    %Storing times
    T(tstep+1) = time;
    
    if filter == 1.0
        if t == 0
            u = F*u;            %Filtering, originally not implemented by default
        elseif t==1
            %u = spdiags(1./sqrt(w'),0,N+1,N+1)*F*(sqrt(W))*u;
            u = F*u;
        elseif t==2
            u = F*u;
        end
    end
    time
    time = time + dt;

    %
    %  Display solution on every 10th time step.
    %
    if ( rem ( tstep, 10 ) == 0 )
%       %ploting analytical
%       plot(x,ua,'k')
%       hold on
      
      for i = 1 : Elements
        plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
          'Markersize', 1, 'LineWidth', 2 );
        hold all
      end
      grid ( 'on' );
      axis ( [ xL, xR, uBott, uUp ] );
      
      drawnow;
      hold off
    end
    pause(0.05)
    
%     %ploting
%     xv = reshape(x,1,[]);
%     uv = reshape(u,1,[]);
%     plot(xv,uv);
%     xlim([xL xR]);
%     ylim([uBott uUp]);
%     drawnow

end

Econ = EEt + dEEt + dfEEt + nfEEt;

plot(T,E)
hold on
plot(T,Ev)
hold on
plot(T,Evf)
hold on
plot(T,Enlf)
hold on
plot(T,E+Ev+Evf+Enlf)
legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
title('Energy vs Time, \nu=0.0, \Omega = whole domain')
%ylim([-E(1)/10 E(1)+E(1)/5])
ylim([-0.1 1.1])
xlim([T(1) T(tstep+1)])
%xlim([T(1) 0.8])
ylabel('Energy')
xlabel('Time')
% 
% a=6;
% liminf = min([EEt(:,a);dEEt(:,a);dfEEt(:,a);nfEEt(:,a)])-EEt(1,a)/10;
% limsup = max([EEt(:,a);dEEt(:,a);dfEEt(:,a);nfEEt(:,a)])+EEt(1,a)/1.5;
% 
% plot(T,EEt(:,a))
% hold on
% plot(T,dEEt(:,a))
% plot(T,dfEEt(:,a))
% plot(T,nfEEt(:,a))
% plot(T,Econ(:,a))
% legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
% title('Energy vs Time, \nu=0.05, \Omega = 6th domain')
% %ylim([liminf limsup])
% %ylim([-0.22 0.32])
% ylim([-0.5 0.2])
% xlim([T(1) T(tstep+1)])
% %xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')

%end Burgers1D subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


