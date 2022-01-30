function [ u ] = Initial_condition ( x , epsilon)
%
% Set this function to evaluate the initial condition of your simulation
%

% u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;
u = -sin(pi*x);

end