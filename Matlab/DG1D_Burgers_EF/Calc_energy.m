%This script computes the Burgers Energy Derived Equation terms inside the
%Runge-Kutta method

Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Derivative of U
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Computing the first derivative (without more nodes)
deru = Dr*u./J;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% u times du in modal (for dissipation fluxes)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
udum = invV*(u.*(Dr*u./J));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% u cube in modal (for non-linear fluxes)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
u3m = invV*(u.^3);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%For viscous fluxes over boundaries
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Over every subdomain
dflux = 2*epsilon.*u([1 N+1],:).*deru([1 N+1],:);
dflux = (dflux(2,:) - dflux(1,:));

%Over every subdomain and mode
dfluxm = zeros(N+1,K);
for i=1:K
    dfm = 2*epsilon.*(udum(:,i).*V([1 N+1],:)');
    dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For non-linear flux over boundaries
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
nflux = 2*(u([1 N+1],:).^3)./3;
nflux = -(nflux(2,:) - nflux(1,:));

% Over every subdomain and mode
nfluxm = zeros(N+1,K);
for i=1:K
    nfm = (2/3)*(u3m(:,i).*V([1 N+1],:)');
    nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For dissipation term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For every subdomain
edis = -(2*epsilon*w*(((deru).^2).*J));

%For every subdomain and mode
edism = -2*epsilon*(invV*(Dr*u./J)).^2.*J(1,:);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% solving energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%RHSE of energy equation integrated over a domain
%rhse = edis; %dflux + nflux + edis;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% dissipation due to viscous term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
rese = rk4a(INTRK)*rese + dt*edis;
EE = EE + rk4b(INTRK)*rese;

% Over every subdomain and mode
resem = rk4a(INTRK)*resem + dt*edism;
EEm = EEm + rk4b(INTRK)*resem;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% viscous fluxes term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
resedf = rk4a(INTRK)*resedf + dt*dflux;
EEdf = EEdf + rk4b(INTRK)*resedf;

% Over every subdomain and mode
resedfm = rk4a(INTRK)*resedfm + dt*dfluxm;
EEdfm = EEdfm + rk4b(INTRK)*resedfm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% non-linear fluxes term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
resenf = rk4a(INTRK)*resenf + dt*nflux;
EEnf = EEnf + rk4b(INTRK)*resenf;

% Over every subdomain and mode
resenfm = rk4a(INTRK)*resenfm + dt*nfluxm;
EEnfm = EEnfm + rk4b(INTRK)*resenfm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





