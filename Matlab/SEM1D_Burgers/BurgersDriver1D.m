% Driver script for solving the 1D burgers equations
%
% PDE:
%           du/dt + du^2/dx = 0  for -1 < x < +1
%         IC:
%           epsilon = 0.1
%           u(x,0) = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0

clear
clc
Globals1D;
%%
%CONTROL PANEL
% Order of polymomials used for approximation
N = 12;

%Final time
FinalTime = 1.2;

%dealiasing??
deal = 1.0;     %1.0 for dealiasing, 0.0 for non-dealiasing

%improved integration??
iint = 0.0;     %1.0 for immproved integration, 0.0 for not
%hay un error con iint=1.0, dan mal los flujos viscosos

%init Filter matrix
%filter = ?
filter = 0;   %1.0 for filtering, 0.0 for non-filtering
t = 2;          % 0 for Non-Unitary, 1 for unitary, 2 for Lobatto Basis.
Nc=1;           %when using Lobatto Bassis you should put Nc=3 for preserve 2 first basis functions (because keep the B.C. of subdomains)
if t==2
    s=16;           % exponential filter degree
else
    s=16;           % exponential filter degree
end

%viscosity constant
%epsilon = 0.01;
nu = 0.0005;
epsilon = nu;

% Generate simple mesh
xL = -1.0;
xR = 1.0;
Elements = 10;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,Elements);
%END CONTROL PANEL

%%

% Initialize solver and construct grid and metric
StartUp1D;

% Filter matrix F
if filter==1.0
    %last parameter (t) means type,
    [F,L] = Filter1D(N,Nc,s,t,V,W,r);
end

% Set initial conditions
%u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;
%uv = exp(-xv.^2/sqrt(0.01));
uv = -sin(pi*xv);
u = vec2matrix(uv,K);
EE = zeros(1,K);
EEdf = zeros(1,K);
EEnf = zeros(1,K);

%Down and Up limits
uBott = min(min(u)) - 0.1;
uUp = max(max(u)) + 0.1;

% Solve Problem
%[u] = Burgers1D ( u, epsilon, xL, xR, FinalTime );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Burgers1D subroutine - Hetshaven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%init time
  time = 0.0;
%
%  Runge-Kutta residual storage  
%
  %resu = zeros(Np,K);
  resu = zeros(Ns,1);
  rese = zeros(1,K);      %result vector for energy disipation due to viscosity
  resedf = zeros(1,K);    %result vector for energy disipation flux over bound.
  resenf = zeros(1,K);    %result vector for energy flux due to non-linear term
%
%  Compute time step size
%
  xmin = min(abs(x(1,:)-x(2,:)));
  CFL = 1.0;
  umax = max(max(abs(u)));
  %dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));     %original
  dt = 0.1*CFL*(xmin/umax);
  dt = 0.0025;
  %dt = CFL*(xmin/umax);
  %dt = 3.1776e-04;
  Nsteps = ceil(FinalTime/dt);
  dt = FinalTime/Nsteps; 

%
%  Matrices
%
  Matrices1D;

%
%  Initial Energy
%
  Init_energy

%
%  Outer time step loop 
%
  figure ( 1 );
  shapestr = { '-o','-x' };
  
for tstep=1:Nsteps
    for INTRK = 1:5
      timelocal = time + rk4c(INTRK)*dt;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Advection step
      
      %[rhsu] = BurgersRHS1D(u,epsilon,xL,xR,timelocal);
      %uv(1) = exp(-(xv(1)-vel*timelocal).^2/sqrt(0.01));
      
      %Advection
        uv(1) = 0;
        u(1) = uv(1);
        uv(Ns) = 0;
        u(end) = uv(Ns);
      %uv(Ns) = exp(-(xv(Ns)-vel*timelocal).^2/sqrt(0.01));
      %rhsu = -1.0*AG*(0.5*uv.^2);

      if deal==1.0
          %Dealiasing is not working well at all
          %umd = invV*u;             %u modal
          umd = [invV*u ; zeros(length(rd)-length(r),Elements)];    %u modal and zero-padding
          ud = Vd*umd;                                              %going to nodal with more nodes
          udv = mat2vec(ud);
          %dfdrd = -0.5*AGd*(udv .^ 2);                                        %computing non-linear term with more nodes
          %dfdrd = vec2matrix(dfdrd,K);
          dfdrd = ud .^ 2;
          umd = invVd*dfdrd;                                        %putting the non-linear result in modal
          umd = umd(1:N+1,:);                                       %extracting the padding part
          ud = V*umd;                                               %going back to original nodal
          uv = mat2vec(ud);

          rhsu = -0.5*AG*(uv.^2);                                   %computing the discretized term
          %rhsu = uv;
      else
          rhsu = -0.5*AG*(uv.^2);
      end
      
        resu = rk4a(INTRK)*resu + dt*rhsu;
        uv = uv+rk4b(INTRK)*resu;

      %end Advection Step
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
      if nu~=0.0
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Diffusion step
      
        ep = 1/((rk4c(INTRK+1)-rk4c(INTRK))*dt*nu);
        uvf=-uv*ep;
        b = zeros(Ns,1);
        DG = SG;
        ap=0; 

        for k=1:K
            for i=1:Np
                b(Np*(k-1)+i-ap)=b(Np*(k-1)+i-ap)-(uvf((Np*(k-1))+i-ap)*w(i)*J(i,k));
            %    Ojo con el signo de aca, pensar bien la vaina
                for j=1:Np
                    DG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) = DG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) + ...
                        ep*M(i,j)*J(i,k);
                end
            end
            ap=ap+1;
        end
      % Boundary conditions
        DG(1,:) = 0.0;
        DG(Ns,:) = 0.0;
        DG(1,1) = 1;
        DG(Ns,Ns) = 1;
        b(1)=0.0;
        b(Ns)=0.0;

      % Solving
        uv = DG\b;
      end
      u = vec2matrix(uv,K);
      
      %end Diffusion Step
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % energy calculations
        
        Calc_energy           
            
      %end energy calculations
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Analitical
    % ua = analitica(x,time+dt,epsilon);
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % energy storage        
    Store_energy
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Filtering     
    if filter==1.0
        uv = filtering(uv,F,K,Np);
    end

%     %Solving diffusion
%     %ep = 1/((rk4c(INTRK+1)-rk4c(INTRK))*dt*nu);
%     uvf=-uv*ep;
%     b = zeros(Ns,1);
%     DG = SG;
%     ap=0; 
%     
%     for k=1:K
%       for i=1:Np
%         b(Np*(k-1)+i-ap)=b(Np*(k-1)+i-ap)-(uvf((Np*(k-1))+i-ap)*w(i)*J(i,k));
%         %    Ojo con el signo de aca, pensar bien la vaina
%         for j=1:Np
%           DG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) = DG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) + ...
%               ep*M(i,j)*J(i,k);
%         end
%       end
%       ap=ap+1;
%     end
%     % Boundary conditions
%     DG(1,:) = 0.0;
%     DG(Ns,:) = 0.0;
%     DG(1,1) = 1;
%     DG(Ns,Ns) = 1;
%     b(1)=0.0;
%     b(Ns)=0.0;
%  
%   % Solving
%     uv = DG\b;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Advancing time 
    time = time + dt;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plotting
    %Getting uv (vector) to v (matrix)
    %u = vec2matrix(uv,K);
    
%     %ploting
    %plot(xv,uv)
    plot(x,u)
    hold on
    %plot(x,ua);
    xlim([xL xR]);
    ylim([-1.1 1.1]);
    grid on
    drawnow
    hold off
%     error = abs(ua-uv);
%     semilogy(xv,error)
%     ylim([1E-10 1E-2])
%     drawnow
    
end

Econ = EEt + dEEt + dfEEt + nfEEt;

hold on
plot(T,E)
hold on
plot(T,Ev)
hold on
plot(T,Evf)
hold on
plot(T,Enlf)
hold on
Ebal = E+Ev+Evf+Enlf;
plot(T,Ebal)
legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
%legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','E1','E_{vd}1','E_{vf}1','E_{nlf}1','E+E_{vd}+E_{vf}+E_{nlf}1','Location','eastoutside')
title('Energy vs Time, \nu=0.0, \Omega = whole domain')
%ylim([-E(1)/10 E(1)+E(1)/5])
ylim([-0.1 1.6])
xlim([T(1) T(tstep+1)])
%xlim([T(1) 0.8])
ylabel('Energy')
xlabel('Time')
% 
% a=5;
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
% title('Energy vs Time, \nu=0.0, \Omega = 5th domain')
% %ylim([liminf limsup])
% %ylim([-0.22 0.32])
% ylim([-0.5 0.2])
% xlim([T(1) T(tstep+1)])
% %xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')

%end Burgers1D subroutine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

