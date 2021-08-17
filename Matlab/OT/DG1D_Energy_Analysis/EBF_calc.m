
Globals1D

%Derivative of U;
deru = Dr*u./J;
%u times du in modal (for dissipation fluxes)
udum = invV*(u.*deru);
%u cube in modal (for non-linear fluxes)
u3m = invV*(u.^3);

%
%%% RHS %%% 
%

%Viscous fluxes
dfluxm = zeros(N+1,K);
for i=1:K
    dfm = 2*epsilon.*(udum(:,i).*V([1 N+1],:)');
    dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
end

%Non-linear fluxes
nfluxm = zeros(N+1,K);
for i=1:K
    nfm = (2/3)*(u3m(:,i).*V([1 N+1],:)');
    nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
end

%Viscous dissipation of Energy
edism = -2*epsilon*(invV*deru).^2.*J(1,:);

%
%%% RK4 %%% 
%

%dissipation due to viscous term
resem = rk4a(INTRK)*resem + dt*edism;
EEm = EEm + rk4b(INTRK)*resem;

%viscous fluxes term
resedfm = rk4a(INTRK)*resedfm + dt*dfluxm;
EEdfm = EEdfm + rk4b(INTRK)*resedfm;

%non-linear fluxes term
resenfm = rk4a(INTRK)*resenfm + dt*nfluxm;
EEnfm = EEnfm + rk4b(INTRK)*resenfm;

