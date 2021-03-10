%This script computes the Burgers Energy Derived Equation terms inside the
%Runge-Kutta method

Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% u^2 in modal (for non-linear fluxes)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
u2m = invV*(u.^2);

%Lobatto basis
u2ml = invL*(u.^2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For advection flux over elements boundaries
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
nflux = a*(u([1 N+1],:).^2);
nflux = -(nflux(2,:) - nflux(1,:));

% Over every subdomain and mode
nfluxm = zeros(N+1,K);
for i=1:K
    nfm = a.*(u2m(:,i).*V([1 N+1],:)');
    nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
end

% Over every subdomain and mode (Lobatto basis)
nfluxml = zeros(N+1,K);
for i=1:K
    nfml = a.*(u2ml(:,i).*L([1 N+1],:)');
    nfluxml(:,i) = -(nfml(:,2) - nfml(:,1));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% solving energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%RHSE of energy equation integrated over a domain
%rhse = edis; %dflux + nflux + edis;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% advection fluxes term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
resenf = rk4a(INTRK)*resenf + dt*nflux;
EEnf = EEnf + rk4b(INTRK)*resenf;

% Over every subdomain and mode
resenfm = rk4a(INTRK)*resenfm + dt*nfluxm;
EEnfm = EEnfm + rk4b(INTRK)*resenfm;

% Over every subdomain and mode (Lobatto basis)
resenfml = rk4a(INTRK)*resenfml + dt*nfluxml;
EEnfml = EEnfml + rk4b(INTRK)*resenfml;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

