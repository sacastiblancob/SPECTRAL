%This script stores the energy calculations (Calc_energy) in the energy
%vectors (Init_energy) for posprocessing comparisons

Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% computing energy before filtering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Et = u.^2;

%computing integral
EEt(tstep+1,:) = w*(Et.*J);

E(tstep+1) = sum(EEt(tstep+1,:));

%
% Dissipation energy total and by mode
%
dEEt(tstep+1,:) = -EE;
dEEtm(:,tstep+1) = reshape(-EEm,1,(N+1)*K);

%
% Dissipation flux energy total and by mode
%
dfEEt(tstep+1,:) = -EEdf;
dfEEtm(:,tstep+1) =  reshape(-EEdfm,1,(N+1)*K);

%
% Nonlinear flux energy total and by mode
%
nfEEt(tstep+1,:) = -EEnf;
nfEEtm(:,tstep+1) = reshape(-EEnfm,1,(N+1)*K);

%
% Total energy in the system viscous, viscous flux, and nonlinear flux
%
Ev(tstep+1) = -1.0*(sum(EE));
Evf(tstep+1) = -1.0*(sum(EEdf));
Enlf(tstep+1) = -1.0*(sum(EEnf));

%    Evtot = sum(EE);
%    Ev(tstep+1) = Ev(tstep) + Evtot;

%Storing times
T(tstep+1) = time;

%
%Energy per mode
%
um = invV*u;
Etm =(invV*u).^2.*J(1,:);
EEtm(:,tstep+1) = reshape(Etm,1,(N+1)*K);

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % storing energy for Hovmoller plots
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% % First term of energy equation u^2 (actual energy of the solution)
% EH(tstep,:) = reshape(u.^2,1,K*(N+1));
% 
% % Second term of energy equation -(2/3)*u^3 (non-linear flux)
% EHnlf(tstep,:) = reshape((2/3)*u.^3,1,K*(N+1));
% 
% % Third term of energy equation 2*nu*u*(du/dx) (viscous flux)
% deru = Dr*u./J;
% EHvf(tstep,:) = reshape(-2*epsilon*u.*deru,1,K*(N+1));
% 
% % Fourth term of energy equation 2*nu*u*(du/dx) (viscous dissipation)
% EHv(tstep,:) = reshape(2*epsilon*(deru.^2),1,K*(N+1));
% 
% % Coordinates
% XCoor(tstep,:) = reshape(x,1,K*(N+1));
% 
% % Times
% Times(tstep,:) = time*ones(1,K*(N+1));
