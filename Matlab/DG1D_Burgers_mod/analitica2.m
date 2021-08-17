%Solver of the analitic solution of the Burgers Ecuation, using
%Gauss-Legendre numerical integration for 11 nodes.
% Analytic solution for Inviscid Burgers' equation for uo = step
% Hesthaven book, page 141
%

function I=analitica2(x,t,a,b,xo)

I = x;
xt = xo + 1.5*t;
I(x<=xt) = a;
I(x>xt) = b;

end