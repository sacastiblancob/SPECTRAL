
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

%
%%% Actual Energy of the solution %%% 
%
invu = invV*u;
Etm =(invu).^2.*J(1,:);

%
% Etm + (-EEm - EEdfm -EEnfm) == Etom  --> By energy conservation
%
%   But, due to numerical round off and truncation errors, energy
%   conservation is not pretty exact, but a tolerance for the difference
%   between LHS and RHS of above equation should be defined.
%

tol = 1;    %Tolerance in percentage of the initial energy.

Etol = abs(sum(Etom))*(tol/100);
LHS = Etm + (-EEm - EEdfm -EEnfm); %Left hand side

Proofd = sum(LHS) < sum(Etom) - Etol;
Proofu = sum(LHS) > sum(Etom) + Etol;
el = 1:K;

EEl = sum(LHS);
elp = el(Proofu);

for i = 1:length(elp)
    Edif = EEl(elp(i)) - EE(elp(i));
    Edd = 0;
    for j = Np:-1:1
        Edd = Edd + Etm(j,elp(i));
        if Edd > Edif
            Edd = Edd - Etm(j,elp(i));
            break
        end
        Etm(j,elp(i)) = 0.0;
    end
    Etm(j,elp(i)) = Etm(j,elp(i)) - (Edif - Edd);
end


%nLHS = Etm + (-EEm - EEdfm -EEnfm);
u = V*(sqrt(abs(Etm./J(1,:))).*sign(invu));


% %u cube in modal (for non-linear fluxes)
% nu3m = invV*(u.^3);



