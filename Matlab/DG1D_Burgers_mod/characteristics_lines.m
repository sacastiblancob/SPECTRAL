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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%CONTROL PANEL
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For analyticial computation of uo=-sin(pi*x)
[etai,wi] = JacobiGQ(0,0,2040);

% Order of polymomials used for approximation
N = 30;

%Final time
FinalTime = 1.2;

%dealiasing??
deal = 0.0;     %1.0 for dealiasing, 0.0 for non-dealiasing

%improved integration??
iint = 0.0;     %1.0 for immproved integration, 0.0 for not
%hay un error con iint=1.0, dan mal los flujos viscosos

%init Filter matrix
%filter = ?
%1.0 for filtering, 0.0 for non-filtering
filter = 0.0;

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
% epsilon = 0.001;
epsilon = 0;

% Generate simple mesh
xL = -1;
xR = 1;
Elements = 8;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,Elements);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%END CONTROL PANEL
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Initialize solver and construct grid and metric
StartUp1D;

% Filter matrix F
if filter==1.0
    %last parameter (t) means type,
    [F,L] = Filter1D(N,Nc,s,t,V,W,r);
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
  CFL = 0.1;
  umax = max(max(abs(u)));
%   dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));     %original
%   dt = 0.1*CFL*(xmin/umax);
%   dt = CFL*(xmin/umax);
%   dt = 0.0005;
  dt = 0.005;
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
  figure ( 1 );
  shapestr = { '-o','-x' };

%For plotting characteristics
xv = reshape(x,1,[]);
chara = zeros(Nsteps,length(xv));
charax = zeros(Nsteps,length(xv));
charat = zeros(Nsteps,length(xv));
  
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

% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         % BurgersRHS1D subroutine
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         %  Define field differences at faces of the K elements.
% %         %
% %           du = zeros(Nfp*Nfaces,K); 
% %           du(:) = u(vmapM) - u(vmapP);
% %         %
% %         %  Impose boundary condition at x=0
% %         %
% %           %uin =-tanh((xL+0.5-time) / (2*epsilon)) + 1.0;   %original
% %           uin = 0;                     %Dirichlet = 0
% % %           uin = uab;                    %Dirichlet for step function
% % %           uin = u((N+1)*Elements);      %Periodic
% %           du(mapI) = 2.0*(u(vmapI)-uin);
% %           %uout=-tanh((xR+0.5-time) / (2*epsilon)) + 1.0;   %original
% %           uout = 0;                    %Dirichlet = 0
% % %           uout = uaa;                   %Dirichlet for step function
% % %           uout = u((N+1)*Elements);     %Periodic
% %           du(mapO) = 2.0*(u(vmapO) - uout);
% %         %
% %         %  Compute Q.
% %         %
% %           q = sqrt ( epsilon ) * ( rx .* ( Dr * u ) - LIFT*(Fscale.*(nx.*du/2.0)));
% %         %
% %         %  Compute jumps DQ at each element face.
% %         %
% %           dq = zeros(Nfaces,K); 
% %           dq(:) = ( q(vmapM) - q(vmapP) ) / 2.0;
% %         %
% %         %  Impose boundary conditions on the jumps.
% %         %  Here, Dirichlet conditions are used.
% %         %
% %           dq(mapI) = 0.0; 
% %           dq(mapO) = 0.0;
% %         %
% %         %  Evaluate DU2, nonlinear flux at faces of the K elements.
% %         %
% %           du2 = zeros(Nfp*Nfaces,K); 
% %           du2(:) = ( u(vmapM).^2 - u(vmapP).^2 ) / 2.0;
% %         % 
% %         %  Impose boundary conditions on the fluxes.
% %         %
% %           du2(mapI) = ( u(vmapI).^2 - uin.^2 ); 
% %           du2(mapO) = ( u(vmapO).^2 - uout.^2 );
% %         %
% %         %  Compute flux
% %         %
% %           maxvel = max ( max ( abs ( u ) ) );
% %         %
% %         %  Penalty scaling -- See Chapter 7.2
% %         %
% % %          tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
% %         %
% %           tau = 0.0;
% %         %
% %         %  Flux term
% %         %
% %           flux = nx .* ( du2 / 2.0 - sqrt ( epsilon ) * dq ) ...
% %             - maxvel / 2.0 .* du - sqrt ( epsilon ) * tau .* du;
% %         %
% %         %  DFDD, local derivatives of field
% %         %
% %           %dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q );       %original
% %           
% %           if deal==1.0
% %               %umd = invV*u;             %u modal
% %               
% %               %u modal and zero-padding
% %               umd = [invV*u ; zeros(length(rd)-length(r),Elements)];
% %               
% %               %going to nodal with more nodes
% %               ud = Vd*umd;
% %               
% %               %computing non-linear term with more nodes
% %               dfdrd = Drd * ( ud .^ 2 / 2);
% %               %dfdrd = ud .^ 2;
% %               
% %               %settling the non-linear result in modal
% %               umd = invVd*dfdrd;
% %               
% %               %substraction the zero padding part
% %               umd = umd(1:N+1,:);
% %               
% %               %going back to original nodal
% %               ud = V*umd;
% % 
% %               %computing the discretized term
% %               dfdr = ud - Dr * (sqrt ( epsilon ) * q );
% %               %dfdr = Dr * ( (ud ./ 2) - sqrt ( epsilon ) * q );  
% %           else
% %               dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q );
% %           end
% %           
% %         %
% %         %  Compute right hand sides of the semi-discrete PDE
% %         %
% %           rhsu = - ( rx .* dfdr - LIFT * ( Fscale .* flux ) );
% %         
% %         %solving velocity  
% %           resu = rk4a(INTRK)*resu + dt*rhsu;
% %           u = u+rk4b(INTRK)*resu;
% %         
% % %           %Ensure continuity at element interfaces
% % %             for k = 2:K
% % %                 u(Np,k-1) = (u(1,k) + u(Np,k-1))/2;
% % %                 u(1,k) = u(Np,k-1);
% % %                 %uf(Np,k-1) = u(Np,k-1);        %uncomment for not filter the boundaries of the elements
% % %                 %uf(1,k) = u(1,k);              %uncomment for not filter the boundaries of the elements
% % %             end
% %           
% %         % end BurgersRHS1D subroutine
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % energy calculations
          Calc_energy             
        %end energy calculations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
%     %computing analytical solution
    ua = analitica3(x,time+dt,epsilon,etai,wi,ai);      %-sin(pi*x)
%     ua = analitica2(x,time+dt,uab,uaa,xao);          %Inviscid Burgers with step
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % energy storage
% Stop condition (sometimes useful)
%     if time > 0.1
%         pause
%     end
    Store_energy
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Filtering    
%     if filter == 1.0
%         u = filtering(u,F);
%         %u = F*u;
%     end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Advancing time  
    time = time + dt;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plotting
    %
    %  Display solution on every 10th time step.
    %
    if ( rem ( tstep, 10 ) == 0 )

      % Plotting analytical and numerical solutions  
%       subplot(1,2,1)
      %ploting analytical
      plot(x,ua,'k')
      hold on
      for ii = 1 : Elements
%         plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
%           'Markersize', 1, 'LineWidth', 2.0, 'Color','b');
        plot ( x(:,ii), u(:,ii), shapestr{1+rem(ii,2)}, ...
          'Markersize', 1, 'LineWidth', 1.5);
        hold all
      end
      grid ( 'on' );
      axis ( [ xL, xR, uBott, uUp ] );
      
      drawnow;
      hold off
      
%      plot ( xv, ustv, 'Markersize', 1, 'LineWidth', 1.5, 'Color','r');
      
%       % Plotting energy by frequency
%       subplot(1,2,2)
%       plot(modes,Etom);
%       xlim([-0.1 N+0.1]);
%       ylim([0 0.2]);
%       legend('1','2','3','4','5','6','7','8','9','10','11','12')
%       drawnow
      
    end
    %pause(0.05)
    
    %ploting
    % xv = reshape(x,1,[]);
    uav = reshape(ua,1,[]);
    chara(tstep,:) = uav;
    charax(tstep,:) = xv;
    charat(tstep,:) = (time-dt)*ones(1,length(xv));

%     subplot(1,2,1)
%     plot(x,u);
%     xlim([xL xR]);
%     ylim([uBott uUp]);
%     drawnow

end
% surf(chara)
% contour3(charax,charat,chara,'k')
% contour3(charax,charat,chara)
% colormap(gray)
contour3(charax,charat,chara,[-0.999,-0.8:0.2:0.8,0.999],'Color',[0.4 0.4 0.4])
view(0,90)
xlabel('X')
ylabel('t')
title(strcat('Characteristics \nu=',string(epsilon)'))
hold on
plot(x,(1/pi)*ones(size(x)),'--k','DisplayName','t_c')
legend('Characteristics','t_c')

% Postprocessing
Econa = EEta + dEEta + dfEEta + nfEEta;
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % Plot of global energy behaviour ANALYTICAL
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(8)

hold on
plot(T,Ea,'LineWidth',1.5)
hold on
plot(T,Eva,'LineWidth',1.5)
hold on
plot(T,Evfa,'LineWidth',1.5)
hold on
plot(T,Enlfa,'LineWidth',1.5)
hold on
Ebala = Ea+Eva+Evfa+Enlfa;
% Eotro = E+Ev;
Eotroa = Ebala;
plot(T,Eotroa,'LineWidth',2)
legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','Location','eastoutside')
legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','E1 a','E_{vd}1 a','E_{vf}1 a','E_{nlf}1 a','E+E_{vd}+E_{vf}+E_{nlf}1 a','Location','eastoutside')
title(strcat('Energy vs Time, \nu=',string(epsilon),', \Omega = whole domain (analytical)'))
ylim([-E(1)/10 E(1)+E(1)/5])
% ylim([-0.9 1.4])
xlim([T(1) T(tstep+1)])
% xlim([T(1) 0.8])
ylabel('Energy (analytical)')
xlabel('Time')
grid on

%Slope of Ebal
EbalEsl = zeros(Nsteps+1,1);
for i=1:Nsteps
    EbalEsl = (Ebal(i+1)-Ebal(i))/dt;
end

% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Plot of global energy behaviour by mode a ANALYTICAL
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(9)

ap = 6;
a = ap;

liminf = min([EEta(:,a);dEEta(:,a);dfEEta(:,a);nfEEta(:,a)])-EEta(1,a)/10;
limsup = max([EEta(:,a);dEEta(:,a);dfEEta(:,a);nfEEta(:,a)])+EEta(1,a)/1.5;

plot(T,EEta(:,a),'LineWidth',1.5)
hold on
plot(T,dEEta(:,a),'LineWidth',1.5)
plot(T,dfEEta(:,a),'LineWidth',1.5)
plot(T,nfEEta(:,a),'LineWidth',1.5)
plot(T,Econa(:,a),'LineWidth',1.5)
legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','Location','eastoutside')
title(strcat('Energy vs Time, \nu =',string(epsilon),', \Omega = ',string(a),'th domain (analytical)'))
%ylim([liminf limsup])
% ylim([-0.22 0.32])
ylim([-0.3 0.4])
% ylim([min([EEt(:,a) ; dEEt(:,a); dfEEt(:,a); nfEEt(:,a); Econ(:,a)])-0.05, max([EEt(:,a) ; dEEt(:,a); dfEEt(:,a); nfEEt(:,a); Econ(:,a)])+0.05])
xlim([T(1) T(tstep+1)])
%xlim([T(1) 0.8])
ylabel('Energy  (analytical)')
xlabel('Time')
grid on

%end Burgers1D subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


