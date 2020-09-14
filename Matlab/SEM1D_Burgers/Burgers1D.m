function [ u ] = Burgers1D ( u, epsilon, xL, xR, FinalTime )

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
      uv(1) = uv(Ns);
      %uv(Ns) = exp(-(xv(Ns)-vel*timelocal).^2/sqrt(0.01));
      %rhsu = -1.0*vel*AG*uv;
      rhsu = -1.0.*uv.*(AG*uv); 
      resu = rk4a(INTRK)*resu + dt*rhsu;
      uv = uv+rk4b(INTRK)*resu;
    end
    time = time + dt;
    
    %ploting
%     xv = reshape(x,1,[]);
%     uv = reshape(u,1,[]);
    plot(xv,uv);
    xlim([xL xR]);
    ylim([-1.1 1.1]);
    drawnow
    
  end

  return
end

