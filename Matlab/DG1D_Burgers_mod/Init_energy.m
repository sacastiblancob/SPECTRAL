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

%3 dimensional tensors for store modal energy
EEtm = zeros(N+1,K,Nsteps+1);   %energy by element and by mode
if iint==0
    dEEtm = zeros(N+1,K,Nsteps+1);  %energy dissipation by element due to viscosity (by mode)
    dfEEtm = zeros(N+1,K,Nsteps+1);  %energy fluxes due to dissipation term (by mode)
    nfEEtm = zeros(N+1,K,Nsteps+1);  %energy fluxes due to non-linear term (by mode)
else
    dEEtm = zeros(Nd+1,K,Nsteps+1);  %energy dissipation by element due to viscosity (by mode)
    dfEEtm = zeros(Nd+1,K,Nsteps+1);  %energy fluxes due to dissipation term (by mode)
    nfEEtm = zeros(Nd+1,K,Nsteps+1);  %energy fluxes due to non-linear term (by mode)
end

%Initial Energy
Eto = u.^2;
Eto = w*(Eto.*J);
EEt(1,:) = Eto;
E(1) = sum(EEt(1,:));
T(1) = 0.0;

%Energía inicial por modo
Etom =(invV*u).^2.*J(1,:);
EEtm(:,:,1) = Etom;

    
% %commented, only for posprocessing purposes
% EM = zeros(Nsteps+1,K);
% for i=1:Nsteps+1
%   EM(i,:) = sum(nfEEtm(:,:,i));
% end




