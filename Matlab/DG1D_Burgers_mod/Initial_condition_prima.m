function [ uprima ] = Initial_condition_prima ( x , epsilon)
%
% Set this function to evaluate the integral from 0 to h of the initial 
% condition of your simulation
%

% u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;
uprima = -pi*cos(pi*x);

end