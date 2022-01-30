% Script to set general parameters
% Sergio Castiblanco
%
% Works under modified 1D-Burgers Equation solver by Prof. Jan Hesthaven
% 

Globals1D;

%%%GENERAL CONFIGURATION


% Order of polymomials used for approximation
N = 15;

% Number of Elements
Elements = 8;

% Boundaries
xL = -1.0;
xR = 1.0;

% Dealiasing? (3/2 rule) 1=yes, 0=no
deal = 0;

% Classical filter? (standard exponential) 1=yes, 0=no
cf = 0;

% SVV approach? 1=yes, 0=no
svv = 1;

% New energy based filter methodology? 1=yes, 0=no
ebf = 0;

% Analytical? 1=yes, 0=no
ana = 1;

%
%%% If Classical exponential filter (cf == 1) then: %%%%%%%%%%%%%%%%%%%%%%%
%
if cf==1
% Order of filter
s = 4;

% Alfa parameter for exponential filter
alp = -log(eps);

% Cutoff for filter (filter will be apply from Nc on)
Nc = 1;
end

%
%%% If SVV approach (svv == 1) then: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if svv==1
% Type of epsilon constant for diffusion term
% eps_svv = 1; standard (1/2)*min(dx)
% eps_svv = 2; -alp*((-1)^2p/(dt*N^(2*p)), being P de order of the SVV term
eps_svv = 1;

% Type of Q operator inside SVV term
% Q_svv = 1;   Standard for Legendre Q = exp(−(k −N)^2/(k −mN)^2) (Pasquetti, 2006)
% Q_svv = 2;   sqrt(1-x^2)
Q_svv = 1;

end

%
%%% If EBF approach (ebf == 1) then: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ebf==1
% Type of filter
% filter_ebf = 1; Exponential
% filter_ebf = 2; Vandeven
filter_ebf = 1;

end



