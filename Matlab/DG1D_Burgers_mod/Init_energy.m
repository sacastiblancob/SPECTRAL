%This script initialize the vectors for store energy results and computes
%the initial energy of the system and by element


Globals1D

%Initialization
E = zeros(Nsteps+1,1);      %Energy of whole domain
Ev = zeros(Nsteps+1,1);     %Energy dissipated by viscosity in whole domain
Evf = zeros(Nsteps+1,1);    %Fluxes of energy due to viscosity in whole domain
Enlf = zeros(Nsteps+1,1);   %Non-linear fluxes of energy in whole domain
dEEt = zeros(Nsteps+1,K);   %energy dissipation by element due to viscosity
dfEEt = zeros(Nsteps+1,K);  %energy flux due to dissipation
nfEEt = zeros(Nsteps+1,K);  %energy flux due to non-linear term
EEt = zeros(Nsteps+1,K);    %energy by element
%Econ = zeros(Nsteps+1,K);   %conservation?
T = zeros(Nsteps+1,1);

%Initial Energy
Eto = u.^2;
EEt(1,:) = w*(Eto.*J);
E(1) = sum(EEt(1,:));
T(1) = 0.0;

%Energía por modo
Etom =(invV*u).^2.*J(1,:);






