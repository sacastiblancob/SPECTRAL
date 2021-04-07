function [x,w] = RWgll(alp,bet,N)

%Computes the roots and Gaussian Quadrature GLL weights for Jacobi
% polynomials if alp = bet in:
%       w(x) = (1-x)^(alp) * (1+x)^(bet)
% if alp =/ beta, there's no human way to compute de GLL weights.

x = JacobiGL(alp,bet,N);
if alp == bet
    w = WGint(alp,N,x);
else
    error('(Jacobi) alpha is not equal beta!!');
end
    
end