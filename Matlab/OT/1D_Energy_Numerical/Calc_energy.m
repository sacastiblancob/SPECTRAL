%This script computes the Burgers Energy Derived Equation terms inside the
%Runge-Kutta method

Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Derivative of U
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if iint==0
    %Computing the first derivative (without more nodes)
    deru = Dr*u./J;
%     derua = Dr*ua./J; %(analytical)
else
   %u modal and zero-padding
    umd = [invV*u ; zeros(length(rd)-length(r),Elements)];
%     umda = [invV*ua ; zeros(length(rd)-length(r),Elements)]; %(analytical)
    
    %to nodal with more modes
    ud = Vd*umd;
%     uda = Vd*umda; %(analytical)
    
    %Computing the first derivative
    deru = Drd*ud./Jd;
%     derua = Drd*uda./Jd; %(analytical)
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% u times du in modal (for dissipation fluxes)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if iint==0
    udum = invV*(u.*(Dr*u./J));
%     uduma = invV*(ua.*(Dr*ua./J)); %(analytical)
% else
%     udum = invVd*(ud.*deru);
%     uduma = invVd*(uda.*deru); %(analytical)
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% u cube in modal (for non-linear fluxes)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if iint==0
    u3m = invV*(u.^3);
%     u3ma = invV*(ua.^3); %(analytical)
% else
%     u3m = invVd*(ud.^3);
%     u3ma = invVd*(uda.^3); %(analytical)
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%For viscous fluxes over boundaries
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Over every subdomain
if iint==0
    dflux = 2*epsilon.*u([1 N+1],:).*deru([1 N+1],:);
    %dflux = 2*epsilon.*u([1 N+1],:);
    dflux = (dflux(2,:) - dflux(1,:));
    
%     % over analytical
%     dfluxa = 2*epsilon.*ua([1 N+1],:).*derua([1 N+1],:);
%     %dflux = 2*epsilon.*u([1 N+1],:);
%     dfluxa = (dfluxa(2,:) - dfluxa(1,:));
else
    dflux = 2*epsilon.*u([1 N+1],:).*deru([1 Nd+1],:);
    dflux = (dflux(2,:) - dflux(1,:));
    
%     % over analytical
%     dfluxa = 2*epsilon.*ua([1 N+1],:).*derua([1 Nd+1],:);
%     dfluxa = (dfluxa(2,:) - dfluxa(1,:));
end
%sum(dflux)

%Over every subdomain and mode
% if iint==0
    dfluxm = zeros(N+1,K);
    for i=1:K
        dfm = 2*epsilon.*(udum(:,i).*V([1 N+1],:)');
        dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
    end
    
%     %over analytical
%     dfluxma = zeros(N+1,K);
%     for i=1:K
%         dfma = 2*epsilon.*(uduma(:,i).*V([1 N+1],:)');
%         dfluxma(:,i) = (dfma(:,2) - dfma(:,1));
%     end
% else
%     dfluxm = zeros(Nd+1,K);
%     for i=1:K
%         dfm = 2*epsilon.*(udum(:,i).*Vd([1 Nd+1],:)');
%         dfluxm(:,i) = (dfm(:,2) - dfm(:,1));
%     end
%     
%     %over analytical
%     dfluxma = zeros(Nd+1,K);
%     for i=1:K
%         dfma = 2*epsilon.*(uduma(:,i).*Vd([1 Nd+1],:)');
%         dfluxma(:,i) = (dfma(:,2) - dfma(:,1));
%     end
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For non-linear flux over boundaries
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
nflux = 2*(u([1 N+1],:).^3)./3;
nflux = -(nflux(2,:) - nflux(1,:));

% nfluxa = 2*(ua([1 N+1],:).^3)./3; %(analytical)
% nfluxa = -(nfluxa(2,:) - nfluxa(1,:)); %(analytical)

% Over every subdomain and mode
% if iint==0
    nfluxm = zeros(N+1,K);
    for i=1:K
        nfm = (2/3)*(u3m(:,i).*V([1 N+1],:)');
        nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
    end
    
%     %over analytical
%     nfluxma = zeros(N+1,K);
%     for i=1:K
%         nfma = (2/3)*(u3ma(:,i).*V([1 N+1],:)');
%         nfluxma(:,i) = -(nfma(:,2) - nfma(:,1));
%     end
% else
%     nfluxm = zeros(Nd+1,K);
%     for i=1:K
%         nfm = (2/3)*(u3m(:,i).*Vd([1 Nd+1],:)');
%         nfluxm(:,i) = -(nfm(:,2) - nfm(:,1));
%     end
%     
%     %over analytical
%     nfluxma = zeros(Nd+1,K);
%     for i=1:K
%         nfma = (2/3)*(u3ma(:,i).*Vd([1 Nd+1],:)');
%         nfluxma(:,i) = -(nfma(:,2) - nfma(:,1));
%     end
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For dissipation term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% For every subdomain
if iint==0
    edis = -(2*epsilon*w*(((deru).^2).*J));
%     edisa = -(2*epsilon*w*(((derua).^2).*J)); %(analytical)
else
    edis = -(2*epsilon*wd*(((deru).^2).*Jd));
%     edisa = -(2*epsilon*wd*(((derua).^2).*Jd)); %(analytical)
end

%For every subdomain and mode
% if iint==0
    edism = -2*epsilon*(invV*(Dr*u./J)).^2.*J(1,:);
%     edisma = -2*epsilon*(invV*(Dr*ua./J)).^2.*J(1,:); %(analytical)
% else
%     edism = -2*epsilon*(invVd*deru).^2.*Jd(1,:);
%     edisma = -2*epsilon*(invVd*derua).^2.*Jd(1,:); %(analytical)
% end

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
% resea = rk4a(INTRK)*resea + dt*edisa; %(analytical)e
% EEa = EEa + rk4b(INTRK)*resea; %(analytical)

% % Over every subdomain and mode
% resem = rk4a(INTRK)*resem + dt*edism;
% EEm = EEm + rk4b(INTRK)*resem;
% resema = rk4a(INTRK)*resema + dt*edisma; %(analytical)
% EEma = EEma + rk4b(INTRK)*resema; %(analytical)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% viscous fluxes term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
resedf = rk4a(INTRK)*resedf + dt*dflux;
EEdf = EEdf + rk4b(INTRK)*resedf;
% resedfa = rk4a(INTRK)*resedfa + dt*dfluxa; %(analytical)
% EEdfa = EEdfa + rk4b(INTRK)*resedfa; %(analytical)

% % Over every subdomain and mode
% resedfm = rk4a(INTRK)*resedfm + dt*dfluxm;
% EEdfm = EEdfm + rk4b(INTRK)*resedfm;
% resedfma = rk4a(INTRK)*resedfma + dt*dfluxma; %(analytical)
% EEdfma = EEdfma + rk4b(INTRK)*resedfma; %(analytical)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% non-linear fluxes term
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Over every subdomain
resenf = rk4a(INTRK)*resenf + dt*nflux;
EEnf = EEnf + rk4b(INTRK)*resenf;
% resenfa = rk4a(INTRK)*resenfa + dt*nfluxa; %(analytical)
% EEnfa = EEnfa + rk4b(INTRK)*resenfa; %(analytical)

% % Over every subdomain and mode
% resenfm = rk4a(INTRK)*resenfm + dt*nfluxm;
% EEnfm = EEnfm + rk4b(INTRK)*resenfm;
% resenfma = rk4a(INTRK)*resenfma + dt*nfluxma; %(analytical)
% EEnfma = EEnfma + rk4b(INTRK)*resenfma; %(analytical)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





