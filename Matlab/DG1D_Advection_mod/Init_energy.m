%This script initialize the vectors for store energy results and computes
%the initial energy of the system and by element


Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Energy Arrays
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%array to store RK4 steps of advection fluxes
EEnf = zeros(1,K);

%Matrices to store results
EEnfm = zeros(N+1,K);     %matrix for store RK4 steps of modal advection fluxes

EEnfml = zeros(N+1,K);     %matrix for store RK4 steps of modal advection fluxes (Lobatto basis)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initialization of Energy Arrays and Matrices
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resu = zeros(Np,K);
resenf = zeros(1,K);    %result vector for energy flux due to advection term

resenfm = zeros(N+1,K);   %result matrix for energy modal advection fluxes

resenfml = zeros(N+1,K);   %result matrix for energy modal advection fluxes (Lobatto basis)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initialization of Energy Arrays and Matrices
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
E = zeros(Nsteps+1,1);      %Energy of whole domain
Enlf = zeros(Nsteps+1,1);   %Advection fluxes of energy in whole domain

nfEEt = zeros(Nsteps+1,K);  %energy flux due advection term
EEt = zeros(Nsteps+1,K);    %energy by element

%Econ = zeros(Nsteps+1,K);   %conservation?
T = zeros(Nsteps+1,1);

%3 dimensional tensors for store modal energy
EEtm = zeros((N+1)*K,Nsteps+1);   %energy by element and by mode

nfEEtm = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to advection term (by mode)

EEtml = zeros((N+1)*K,Nsteps+1);   %energy by element and by mode (Lobatto basis)

nfEEtml = zeros((N+1)*K,Nsteps+1);  %energy fluxes due to advection term (by mode) (Lobatto basis)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initial Energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Eto = u.^2;
Eto = w*(Eto.*J);
EEt(1,:) = Eto;
E(1) = sum(EEt(1,:));
T(1) = 0.0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Initial Energy per mode
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Etom =(invV*u).^2.*J(1,:);
EEtm(:,1) = reshape(Etom,1,(N+1)*K);

%Lobatto basis

uml = invL*(u.^2);
Etoml = zeros(N+1,K);
for i = 1:K
   lode = L.*uml(:,i)';
   inte = sum(w'.*(lode),1);
   Etoml(:,i) = inte'.*J(:,i);
end
EEtml(:,1) = reshape(Etoml,1,(N+1)*K);
Etol = sum(Etoml);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Matrices for Hovmoller plots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% First term of energy equation u^2
EH = zeros(Nsteps,K*(N+1));

% Second term of energy equation -a*u^2
EHnlf = zeros(Nsteps,K*(N+1));

% Matrix of coordinates
XCoor = zeros(Nsteps,K*(N+1));

% Matrix of times
Times = zeros(Nsteps,K*(N+1));

%%%%%%

modes = (0:N)';

%end
