function [ uinte ] = Initial_condition_int ( x , epsilon)
%
% Set this function to evaluate the integral from some constant value to h 
% of the initial condition of your simulation.
%
% Usually the constant part of the result of the integral cancels out into
% the evaluation of the analityc solution expression.
%
% See Cole 1951 - On a quasi-linear parabolic equation occurring in aerodynamics
% Page 231 of (225-236)
%

% uinte = x - 2 * epsilon * log( cosh ( ( 0.5*x + 0.25 ) / ( epsilon ) )) + ...
%     2*epsilon*log( cosh ( 0.25 / epsilon ));
uinte = cos(pi*x);

end