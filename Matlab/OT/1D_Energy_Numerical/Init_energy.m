%This script initialize the vectors for store energy results and computes
%the initial energy of the system and by element


Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Energy Arrays
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%array to store RK4 steps of viscous dissipation
EE = zeros(1,K);
EEa = zeros(1,K);

%array to store RK4 steps of viscous fluxes
EEdf = zeros(1,K);
EEdfa = zeros(1,K);

%array to store RK4 steps of non-linear fluxes
EEnf = zeros(1,K);
EEnfa = zeros(1,K);

%Matrices to store results
% if iint==0
    EEm = zeros(N+1,K);     %matrix for store RK4 steps of modal viscous dissipation
    EEdfm = zeros(N+1,K);     %matrix for store RK4 steps of modal viscous fluxes
    EEnfm = zeros(N+1,K);     %matrix for store RK4 steps of modal non-linear fluxes
    
    %analytical
    EEma = zeros(N+1,K);     %matrix for store RK4 steps of modal viscous dissipation
    EEdfma = zeros(N+1,K);     %matrix for store RK4 steps of modal viscous fluxes
    EEnfma = zeros(N+1,K);     %matrix for store RK4 steps of modal non-linear fluxes
% else
%     EEm = zeros(Nd+1,K);     %matrix for store RK4 steps of modal viscous dissipation
%     EEdfm = zeros(Nd+1,K);     %matrix for store RK4 steps of modal viscous fluxes
%     EEnfm = zeros(Nd+1,K);     %matrix for store RK4 steps of modal non-linear fluxes
%     
%     %analytical
%     EEma = zeros(Nd+1,K);     %matrix for store RK4 steps of modal viscous dissipation
%     EEdfma = zeros(Nd+1,K);     %matrix for store RK4 steps of modal viscous fluxes
%     EEnfma = zeros(Nd+1,K);     %matrix for store RK4 steps of modal non-linear fluxes
% end
    
%Down and Up limits
uBott = min(min(u)) - 0.1;
uUp = max(max(u)) + 0.1;

% Solve Problem
%[u] = Burgers1D ( u, epsilon, xL, xR, FinalTime );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Burgers1D subroutine - Hetshaven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%init time
time = 0.0;
%
%  Runge-Kutta residual storage  
%
resu = zeros(Np,K);
rese = zeros(1,K);      %result vector for energy dissipation due to viscosity
resedf = zeros(1,K);    %result vector for energy dissipation flux over bound.
resenf = zeros(1,K);    %result vector for energy flux due to non-linear term

%analytical
resua = zeros(Np,K);
resea = zeros(1,K);      %result vector for energy dissipation due to viscosity
resedfa = zeros(1,K);    %result vector for energy dissipation flux over bound.
resenfa = zeros(1,K);    %result vector for energy flux due to non-linear term


% if iint==0
  resem = zeros(N+1,K);   %result matrix for energy modal disipation due to viscosity
  resedfm = zeros(N+1,K);   %result matrix for energy modal disipation fluxes
  resenfm = zeros(N+1,K);   %result matrix for energy modal non-linear fluxes
  
  %analytical
  resema = zeros(N+1,K);   %result matrix for energy modal disipation due to viscosity
  resedfma = zeros(N+1,K);   %result matrix for energy modal disipation fluxes
  resenfma = zeros(N+1,K);   %result matrix for energy modal non-linear fluxes
% else
%   resem = zeros(Nd+1,K);   %result matrix for energy modal disipation due to viscosity
%   resedfm = zeros(Nd+1,K);   %result matrix for energy modal disipation fluxes
%   resenfm = zeros(Nd+1,K);   %result matrix for energy modal non-linear fluxes
%   
%   %analytical
%   resema = zeros(Nd+1,K);   %result matrix for energy modal disipation due to viscosity
%   resedfma = zeros(Nd+1,K);   %result matrix for energy modal disipation fluxes
%   resenfma = zeros(Nd+1,K);   %result matrix for energy modal non-linear fluxes
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initialization of Energy Arrays and Matrices
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
E = zeros(Nsteps+1,1);      %Energy of whole domain
Ea = zeros(Nsteps+1,1);      %Energy of whole domain fo (analytical)
Ev = zeros(Nsteps+1,1);     %Energy dissipated by viscosity in whole domain
Eva = zeros(Nsteps+1,1);     %Energy dissipated by viscosity in whole domain (analytical)
Evf = zeros(Nsteps+1,1);    %Fluxes of energy due to viscosity in whole domain
Evfa = zeros(Nsteps+1,1);    %Fluxes of energy due to viscosity in whole domain (analytical)
Enlf = zeros(Nsteps+1,1);   %Non-linear fluxes of energy in whole domain
Enlfa = zeros(Nsteps+1,1);   %Non-linear fluxes of energy in whole domain (analytical)
dEEt = zeros(Nsteps+1,K);   %energy dissipation by element due to viscosity
dEEta = zeros(Nsteps+1,K);   %energy dissipation by element due to viscosity (analytical)
dfEEt = zeros(Nsteps+1,K);  %energy flux due to dissipation
dfEEta = zeros(Nsteps+1,K);  %energy flux due to dissipation (analytical)
nfEEt = zeros(Nsteps+1,K);  %energy flux due to non-linear term
nfEEta = zeros(Nsteps+1,K);  %energy flux due to non-linear term (analytical)
EEt = zeros(Nsteps+1,K);    %energy by element
EEta = zeros(Nsteps+1,K);    %energy by element (analytical)
%Econ = zeros(Nsteps+1,K);   %conservation?
T = zeros(Nsteps+1,1);

%3 dimensional tensors for store modal energy
EEtm = zeros((N+1)*K,Nsteps+1);   %energy by element and by mode
EEtma = zeros((N+1)*K,Nsteps+1);   %energy by element and by mode (analytical)
% if iint==0
    dEEtm = zeros((N+1)*K,Nsteps+1);  %energy dissipation by element due to viscosity (by mode)
    dfEEtm = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to dissipation term (by mode)
    nfEEtm = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to non-linear term (by mode)
    
    dEEtma = zeros((N+1)*K,Nsteps+1);  %energy dissipation by element due to viscosity (by mode) (analytical)
    dfEEtma = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to dissipation term (by mode) (analytical)
    nfEEtma = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to non-linear term (by mode) (analytical)
% else
%     dEEtm = zeros((Nd+1)*K,Nsteps+1);  %energy dissipation by element due to viscosity (by mode)
%     dfEEtm = zeros((Nd+1)*K,Nsteps+1);  %energy fluxes due to dissipation term (by mode)
%     nfEEtm = zeros((Nd+1)*K,Nsteps+1);  %energy fluxes due to non-linear term (by mode)
%     
%     dEEtma = zeros((N+1)*K,Nsteps+1);  %energy dissipation by element due to viscosity (by mode) (analytical)
%     dfEEtma = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to dissipation term (by mode) (analytical)
%     nfEEtma = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to non-linear term (by mode) (analytical)
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initial Energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Eto = u.^2;
Eto = w*(Eto.*J);
EEt(1,:) = Eto;
E(1) = sum(EEt(1,:));
T(1) = 0.0;

Ea(1) = sum(EEt(1,:)); %analytical
EEta = Eto;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initial Energy per mode
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Etom =(invV*u).^2.*J(1,:);
EEtm(:,1) = reshape(Etom,1,(N+1)*K);
EEtma(:,1) = reshape(Etom,1,(N+1)*K);

    
% %commented, only for posprocessing purposes
% EM = zeros(Nsteps+1,K);
% for i=1:Nsteps+1
%   EM(i,:) = sum(nfEEtm(:,:,i));
% end

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %Matrices for Hovmoller plots
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% % First term of energy equation u^2
% EH = zeros(Nsteps,K*(N+1));
% 
% % Second term of energy equation -(2/3)*u^3
% EHnlf = zeros(Nsteps,K*(N+1));
% 
% % Third term of energy equation 2*nu*u*(du/dx)
% EHvf = zeros(Nsteps,K*(N+1));
% 
% % Fourth term of energy equation 2*nu*u*(du/dx)
% EHv = zeros(Nsteps,K*(N+1));
% 
% % Matrix of coordinates
% XCoor = zeros(Nsteps,K*(N+1));
% 
% % Matrix of times
% Times = zeros(Nsteps,K*(N+1));
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % end
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










