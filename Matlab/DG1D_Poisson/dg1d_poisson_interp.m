function value = dg1d_poisson_interp ( x, i, nel, order, c )

%*****************************************************************************80
%
%% DG1D_POISSON_INTERP evaluates a DG interpolant at a point in a subinterval.
%
%  Modified:
%
%    16 September 2018
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X, the location of the point.
%
%    Input, integer I, the subinterval containing X.
%
%    Input, integer NE, the number of subintervals.
%
%    Input, integer ORDER, the number of monomials used.
%
%    Input, real C(ORDER*NE), the interpolant coefficients.
%
%    Output, real VALUE, the value of the interpolant at X.
%
  value = 0.0;
  for k = 1 : order
    degree = k - 1 ;
    value = value + c(k+(i-1)*order) * dg1d_poisson_monomial ( x, i, nel, degree );
  end

  return
end
