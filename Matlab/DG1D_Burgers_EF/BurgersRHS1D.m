function [ rhsu ] = BurgersRHS1D ( u, epsilon, xL, xR, time, etai, wi, ai )

%*****************************************************************************80
%
%% BURGERSRHS1D evaluates RHS flux in 1D viscous Burgers equation.
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
%
%  Define field differences at faces of the K elements.
%
  du = zeros(Nfp*Nfaces,K); 
  du(:) = u(vmapM) - u(vmapP);
%
%  Impose boundary condition at x=0
%
%   uin =-tanh((xL+0.5-time) / (2*epsilon)) + 1.0;
  uin = 0;
%   uin = u((N+1)*Elements);
%   if(time == 0)
%       uin = u(1);
%   else
%       uin = analitica(xL,time,epsilon,etai,wi,ai);
%   end
  du(mapI) = 2.0*(u(vmapI)-uin);
%   uout=-tanh((xR+0.5-time) / (2*epsilon)) + 1.0;
  uout = 0;
%   uout = u((N+1)*Elements);
%   if(time == 0)
%       uout = u(end);
%   else
%       uout = analitica(xR,time,epsilon,etai,wi,ai);
%   end
  du(mapO) = 2.0*(u(vmapO) - uout);
%
%  Compute Q.
%
  if (svv == 1)
      epsilon_new_sq = (epsilon + Qsvv*Esvv).^0.5;
      qaux = ( rx .* ( Dr * u ) - LIFT*(Fscale.*(nx.*du/2.0)));
      qauxm = invV*qaux;
      qauxm = epsilon_new_sq.*qauxm;
      q = V*qauxm;
  else
      q = sqrt ( epsilon ) * ( rx .* ( Dr * u ) - LIFT*(Fscale.*(nx.*du/2.0)));
  end
%
%  Compute jumps DQ at each element face.
%
  if (svv == 1)
      qauxm_2 = epsilon_new_sq.*qauxm;
      sqe_q = V*qauxm_2;
      dqe_dq = zeros(Nfaces,K);
      dqe_dq(:) = ( sqe_q(vmapM) - sqe_q(vmapP) ) / 2.0;
      
      dqe_dq(mapI) = 0.0; 
      dqe_dq(mapO) = 0.0;
      
      eps_n_nodal = V*epsilon_new_sq;
      d_eps_n = zeros(Nfaces,K);
      d_eps_n(:) = ( eps_n_nodal(vmapM) - eps_n_nodal(vmapP) ) / 2.0;
  else
      dq = zeros(Nfaces,K); 
      dq(:) = ( q(vmapM) - q(vmapP) ) / 2.0;
%
%  Impose boundary conditions on the jumps.
%  Here, Dirichlet conditions are used.
%
      dq(mapI) = 0.0; 
      dq(mapO) = 0.0;
  end
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
  if (svv == 1)
    flux = nx .* ( du2 / 2.0 - dqe_dq ) ...
        - maxvel / 2.0 .* du - abs(d_eps_n).* tau .* du;
%     flux = nx .* ( du2 / 2.0 - sqrt ( epsilon ) * dq ) ...
%         - maxvel / 2.0 .* du - sqrt ( epsilon ) * tau .* du;
  else
    flux = nx .* ( du2 / 2.0 - sqrt ( epsilon ) * dq ) ...
        - maxvel / 2.0 .* du - sqrt ( epsilon ) * tau .* du;
  end
%
%  DFDD, local derivatives of field
%
%   dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q );

  if (dealcond == 1)
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
      if (svv == 1)
          sqrte_q = V*(qauxm.*epsilon_new_sq);
          dfdr = Dr * ( u .^ 2 / 2 - sqrte_q );
      else
          dfdr = Dr * ( u .^ 2 / 2 - sqrt ( epsilon ) * q );
      end
  end
%
%  Compute right hand sides of the semi-discrete PDE
%
  rhsu = - ( rx .* dfdr - LIFT * ( Fscale .* flux ) );
  return
end

