% function [ u ] = Burgers1D ( u, epsilon, xL, xR, FinalTime )

%*****************************************************************************80
%
%% BURGERS1D integrates the 1D Burgers equation until the final time.
%
%  Discussion:
%
%    An initial condition u in the domain [xL,xR] is given.
%    sqrt(epsilon) is the coefficient of viscosity
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
  Globals1D;

  time = 0.0;
%
%  Compute time step size
%
  xmin = min(abs(x(1,:)-x(2,:)));
  CFL = 0.25;
  umax = max(max(abs(u)));
  dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));
  dt = 5*dt;
  Nsteps = ceil(FinalTime/dt);
  dt = FinalTime/Nsteps;
%
% Initializing Runge-Kutta residual storage and Energy
%
  Init_energy
  modes = (0:N)';
%
% For analytical
%
  if epsilon > 0.001
    ai = -3;
  else
    ai = -3;
  end

  [etai,wi] = JacobiGQ(0,0,2048);

%
%  If classical filter
%
if (cf==1)
    [F,L] = Filter1D(N,Nc,s,V,W,alp);
end

%
%  If SVV
%
if (svv == 1)
    %Epsilon of SVV (O(L/2N))
%     Esvv = 0.5*min(x(2:N+1,:) - x(1:N,:));
%     Esvv = 0.5*(1/(N+1))*(VX(2:K+1)-VX(1:K));
%     Esvv = (1/(N^2))*(VX(2:K+1)-VX(1:K));
    if eps_svv == 1
        Esvv = 0.5*min(x(2:N+1,:) - x(1:N,:));
    elseif eps_svv == 2
%         Esvv = -log(eps())*(1/(dt*(N^2)))*ones(1,K);
%         Esvv = (dt/(N^2))*ones(1,K);
%         Esvv = -log(eps())*(1/(dt*((K*N)^2)))*J(1,:);
        Esvv = -log(eps())*(1/(dt*((K*N)^2)))*min(x(2:N+1,:) - x(1:N,:));
    end
    
    % Q's of SVV
    if Q_svv == 1
        Qsvv = zeros(N+1,1);
        for j=1:N+1
            if j>sqrt(N+1)
                Qsvv(j) = exp(-1*((N+1-j)/(sqrt(N+1)-j))^2);
            end
        end
    elseif Q_svv == 2
        Qsvv = zeros(N+1,1);
        for j=1:N+1
            xaux = 1-j/(N+1);
            Qsvv(j) = sqrt(1 - xaux^2);
        end
    end
end

%
%  If Energy Based Methodology
%
if (ebf == 1)
    %Filterdiags
    Fdiags = ones(N+1,K);
    
    Nc0 = 3*N/4;
    Nc0 = round((Nc0)/2)*2;
    s0 = N/2;
    
    knaf = ones(1,K);
    
    %Nc's
    Ncs = (N+1)*ones(1,K);
%     Ncs = Nc0*ones(1,K);

    Eslopew0 = 0;
    Eslope0 = zeros(1,K);
end

%
% Error vectors
%
  Error = zeros(size(u));
  Err1_t = zeros(Nsteps+1,1);
  Err2_t = zeros(Nsteps+1,1);
  Errinf_t = zeros(Nsteps+1,1);
  Errord_t = log10(dt^4)*ones(Nsteps+1,1);

% U limits for plotting
  ulup = max(max(u));
  ulbo = min(min(u));

%
%  Outer time step loop 
%
  for tstep=1:Nsteps

    for INTRK = 1:5
      timelocal = time + rk4c(INTRK)*dt;
      [rhsu] = BurgersRHS1D(u,epsilon,xL,xR,timelocal,etai,wi,ai);
      resu = rk4a(INTRK)*resu + dt*rhsu;
      u = u+rk4b(INTRK)*resu;
      
      %Computing Energy
      Calc_energy
    end
    
    % Storing Energy
    Store_energy
    
    % Applying Energy Based Filter Methodology, if so
    if ebf == 1
        EB_filter
    end
    
    % Applying Classical Filter, if so
    if cf == 1
        u = F*u;
    end
    
    %Computing analytic
    if ana==1
      ua = analitica(x,time+dt,epsilon,etai,wi,ai); %Analytical solution
    end
    time = time + dt;
    
    %ploting
    if ana == 1
        xv = reshape(x,1,[]);
        uv = reshape(u,1,[]);
        uva = reshape(ua,1,[]);
        
        Error = abs(ua - u);
        Err1_t(tstep+1) = norm(reshape(abs(ua(2:N,:)-u(2:N,:)),1,[]),1);
        Err2_t(tstep+1) = norm(reshape(abs(ua(2:N,:)-u(2:N,:)),1,[]),2);
        Errinf_t(tstep+1) = norm(reshape(abs(ua(2:N,:)-u(2:N,:)),1,[]),'Inf');
        Errord_t(tstep+1) = mean(mean(log10(abs(ua - u))));
        Errorv = reshape(Error,1,[]);
        
        figure(1)
        
        subplot(1,2,1)
        plot(xv,uv);
        hold on
        plot(xv,uva)
        xlim([xL xR]);
        ylim([ulbo-0.1 ulup+0.1]);
        title('Solution u, ua')
        xlabel('X')
        ylabel('U')
        hold off
        
        subplot(1,2,2)
        semilogy(xv,Errorv);
        xlim([xL xR]);
        ylim([1E-16 1E-1]);
        title('Error Structure |ua - u|')
        xlabel('X')
        ylabel('log(Error)')
        
        drawnow
    else
        figure(1)
        xv = reshape(x,1,[]);
        uv = reshape(u,1,[]);
        plot(xv,uv);
        hold on
        xlim([xL xR]);
        ylim([ulbo-0.1 ulup+0.1]);
        title('Solution u')
        xlabel('X')
        ylabel('U')
        
        drawnow
        hold off
    end
        
  end

%   return
% end

