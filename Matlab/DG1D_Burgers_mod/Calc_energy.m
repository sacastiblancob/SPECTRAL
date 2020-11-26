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

%u times du in modal (for dissipation fluxes)
if iint==0
    udum = invV*(u.*deru);
else
    udum = invVd*(ud.*deru);
end

%u cube in modal (for non-linear fluxes)
if iint==0
    u3m = invV*(u.^3);
else
    u3m = invVd*(ud.^3);
end

%For viscous fluxes over boundaries
if iint==0
    dflux = 2*epsilon.*u([1 N+1],:).*deru([1 N+1],:);
    %dflux = 2*epsilon.*u([1 N+1],:);
    dflux = (dflux(2,:) - dflux(1,:));
else
    dflux = 2*epsilon.*u([1 N+1],:).*deru([1 Nd+1],:);
    dflux = (dflux(2,:) - dflux(1,:));
end
sum(dflux)

if iint==0
    dfluxm = zeros(N+1,K);
    for i=1:K
        dfm = 2*epsilon.*(udum(:,i).*V([1 N+1],:)');
        dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
    end
else
    dfluxm = zeros(Nd+1,K);
    for i=1:K
        dfm = 2*epsilon.*(udum(:,i).*Vd([1 Nd+1],:)');
        dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
    end
end

%For non-linear flux over boundaries
nflux = 2*(u([1 N+1],:).^3)./3;
nflux = -(nflux(2,:) - nflux(1,:));

if iint==0
    nfluxm = zeros(N+1,K);
    for i=1:K
        nfm = (2/3)*(u3m(:,i).*V([1 N+1],:)');
        nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
    end
else
    nfluxm = zeros(Nd+1,K);
    for i=1:K
        nfm = (2/3)*(u3m(:,i).*Vd([1 Nd+1],:)');
        nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
    end
end

%For dissipation term
if iint==0
    edis = -(2*epsilon*w*(((deru).^2).*J));
    edism = -2*epsilon*(invV*deru).^2.*J(1,:);
else
    edis = -(2*epsilon*wd*(((deru).^2).*Jd));
    edism = -2*epsilon*(invVd*deru).^2.*Jd(1,:);
end

%solving energy
%RHSE of energy equation integrated over a domain
%rhse = edis; %dflux + nflux + edis;

%dissipation due to viscous term
rese = rk4a(INTRK)*rese + dt*edis;
EE = EE + rk4b(INTRK)*rese;

resem = rk4a(INTRK)*resem + dt*edism;
EEm = EEm + rk4b(INTRK)*resem;

%viscous fluxes term
resedf = rk4a(INTRK)*resedf + dt*dflux;
EEdf = EEdf + rk4b(INTRK)*resedf;

resedfm = rk4a(INTRK)*resedfm + dt*dfluxm;
EEdfm = EEdfm + rk4b(INTRK)*resedfm;

%non-linear fluxes term
resenf = rk4a(INTRK)*resenf + dt*nflux;
EEnf = EEnf + rk4b(INTRK)*resenf;

resenfm = rk4a(INTRK)*resenfm + dt*nfluxm;
EEnfm = EEnfm + rk4b(INTRK)*resenfm;



