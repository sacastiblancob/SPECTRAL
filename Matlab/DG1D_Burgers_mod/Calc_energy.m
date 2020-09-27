%This script computes the Burgers Energy Derived Equation terms inside the
%Runge-Kutta method

Globals1D

%Derivative of U
if iint==0
    deru = Dr*u./J;        %without more nodes
else
    umd = [invV*u ; zeros(length(rd)-length(r),Elements)];    %u modal and zero-padding
    ud = Vd*umd;                                              %going to nodal with more nodes
    deru = Drd*ud./Jd;
end

%u and its derivative in modal
derum = invV*deru;      %derivative of u in modal
um = invV*u;            %u in modal

%For viscous fluxes over boundaries
if iint==0
    dflux = 2*epsilon.*u([1 N+1],:).*deru([1 N+1],:);
    dflux = (dflux(2,:) - dflux(1,:));
else
    dflux = 2*epsilon.*u([1 N+1],:).*deru([1 Nd+1],:);
    dflux = (dflux(2,:) - dflux(1,:));
end

dfluxm = zeros(N+1,K);
for i=1:K
    %dfm = 2*epsilon.*(um(:,i).*V([1 N+1],:)').*(derum(:,i).*V([1 N+1],:)');
    dfm = 2*epsilon.*(um(:,i).*derum(:,i)).*V([1 N+1],:)';
    dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
end

%For non-linear flux over boundaries
nflux = 2*(u([1 N+1],:).^3)./3;
nflux = -(nflux(2,:) - nflux(1,:));

%For dissipation term
if iint==0
    edis = -(2*epsilon*w*(((deru).^2).*J));
    edism = 2*epsilon*(invV*deru).^2.*J(1,:);
else
    edis = -(2*epsilon*wd*(((deru).^2).*Jd));
    edism = 2*epsilon*(invVd*deru).^2.*Jd(1,:);
end

%solving energy
%RHSE of energy equation integrated over a domain
%rhse = edis; %dflux + nflux + edis;

%dissipation due to viscous term
rese = rk4a(INTRK)*rese + dt*edis;
EE = EE + rk4b(INTRK)*rese;

%viscous fluxes term
resedf = rk4a(INTRK)*resedf + dt*dflux;
EEdf = EEdf + rk4b(INTRK)*resedf;

%non-linear fluxes term
resenf = rk4a(INTRK)*resenf + dt*nflux;
EEnf = EEnf + rk4b(INTRK)*resenf;  