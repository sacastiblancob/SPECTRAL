function [ dP ] = GradJacobiP ( r, alpha, beta, N )

%*****************************************************************************80
%
%% GRADJACOBIP evaluates the derivative of the Jacobi polynomial.
%
%  Discussion:
%
%    The Jacobi polynomial is of type (alpha,beta)>-1.
%	 Value are computed at points r for order N and returned in dP[1:length(r))]        
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
%    17 September 2018
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
  dP = zeros(length(r), 1);

  if(N == 0)
    dP(:,:) = 0.0;
  else
    dP = sqrt(N*(N+alpha+beta+1)) * JacobiP(r(:),alpha+1,beta+1, N-1);
  end

  return
end
