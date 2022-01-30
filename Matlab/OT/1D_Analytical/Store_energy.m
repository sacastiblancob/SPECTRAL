%This script stores the energy calculations (Calc_energy) in the energy
%vectors (Init_energy) for posprocessing comparisons

Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% computing energy before filtering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if iint == 0
%     Et = u.^2;
    Eta = ua.^2; %(analytical)

    %computing integral with default nodes
%     EEt(tstep+1,:) = w*(Et.*J);
    EEta(tstep+1,:) = w*(Eta.*J); %(analytical)
else
    %u modal and zero-padding
%     umd = [invV*u ; zeros(length(rd)-length(r),Elements)];
    umda = [invV*ua ; zeros(length(rd)-length(r),Elements)]; %(analytical)
    
    %to nodal with more modes
%     ud = Vd*umd;
    uda = Vd*umda; %(analytical)
    
    %squaring
%     Et = ud.^2;
    Eta = uda.^2; %(analytical)

    %computing integral with more nodes
%     EEt(tstep+1,:) = wd*(Et.*Jd);
    EEta(tstep+1,:) = wd*(Eta.*Jd); %(analytical)
%     umd = umd*0.0;
%     ud = ud*0.0;
    umda = umda*0.0; %(analytical)
    uda = uda*0.0; %(analytical)
end

% E(tstep+1) = sum(EEt(tstep+1,:));
Ea(tstep+1) = sum(EEta(tstep+1,:));  %(analytical)

%     %Computing dissipate energy due to viscosity
%     %duv = (1./J).*Dr*u;
%     dus = (Dr*u./J).^2;
%     %function of N, Elements, CFL, xL and xR
%     %Evele = 2*0.008059650327027*epsilon*w*(dus.*J);
%     Evele = dt*epsilon*w*(dus.*J);
%     %Evele = epsilon*w*(dus.*J);
%
% Dissipation energy total and by mode
%
% dEEt(tstep+1,:) = -EE;
% dEEtm(:,tstep+1) = reshape(-EEm,1,(N+1)*K);

%analytical
dEEta(tstep+1,:) = -EEa;
% dEEtma(:,tstep+1) = reshape(-EEma,1,(N+1)*K);

%
% Dissipation flux energy total and by mode
%
% dfEEt(tstep+1,:) = -EEdf;
% dfEEtm(:,tstep+1) =  reshape(-EEdfm,1,(N+1)*K);

%analytical
dfEEta(tstep+1,:) = -EEdfa;
% dfEEtma(:,tstep+1) =  reshape(-EEdfma,1,(N+1)*K);

%
% Nonlinear flux energy total and by mode
%
% nfEEt(tstep+1,:) = -EEnf;
% nfEEtm(:,tstep+1) = reshape(-EEnfm,1,(N+1)*K);

%analytical
nfEEta(tstep+1,:) = -EEnfa;
% nfEEtma(:,tstep+1) = reshape(-EEnfma,1,(N+1)*K);

%
% Total energy in the system viscous, viscous flux, and nonlinear flux
%
% Ev(tstep+1) = -1.0*(sum(EE));
% Evf(tstep+1) = -1.0*(sum(EEdf));
% Enlf(tstep+1) = -1.0*(sum(EEnf));

%analytical
Eva(tstep+1) = -1.0*(sum(EEa));
Evfa(tstep+1) = -1.0*(sum(EEdfa));
Enlfa(tstep+1) = -1.0*(sum(EEnfa));

%    Evtot = sum(EE);
%    Ev(tstep+1) = Ev(tstep) + Evtot;

%Storing times
T(tstep+1) = time;

% %
% %Energy per mode
% %
% um = invV*u;
% Etm =(invV*u).^2.*J(1,:);
% EEtm(:,tstep+1) = reshape(Etm,1,(N+1)*K);
% 
% %analytical
% uma = invV*ua;
% Etma =(invV*ua).^2.*J(1,:);
% EEtma(:,tstep+1) = reshape(Etma,1,(N+1)*K);

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Energy Based Filter, first try
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% tol = 1;    %Tolerance in percentage of the initial energy.
% 
% Etol = abs(sum(Etom))*(tol/100);
% LHS = Etm + (-EEm - EEdfm -EEnfm); %Left hand side
% 
% Proofd = sum(LHS) < sum(Etom) - Etol;
% Proofu = sum(LHS) > sum(Etom) + Etol;
% el = 1:K;
% 
% EEl = sum(LHS);
% elp = el(Proofu);
% 
% for i = 1:length(elp)
% %    Edif = (EEl(elp(i)) - EE(elp(i)))^0.5;
%     Edif = (EEl(elp(i)) - EE(elp(i)));
%     Edd = 0;
%     for j = Np:-1:1
%         Edd = Edd + Etm(j,elp(i));
%         if Edd > Edif
%             Edd = Edd - Etm(j,elp(i));
%             break
%         end
%         Etm(j,elp(i)) = 0.0;
%     end
%     Etm(j,elp(i)) = Etm(j,elp(i)) - (Edif - Edd);
% end
% 
% 
% %nLHS = Etm + (-EEm - EEdfm -EEnfm);
% invu = invV*u;
% u = V*(sqrt(abs(Etm./J(1,:))).*sign(invu));
% 
% for k = 2:K
%     u(Np,k-1) = (u(1,k)/2) + (u(Np,k-1)/2);
%     %uf(Np,k-1) = u(Np,k-1);        %uncomment for not filter the boundaries of the elements
%     %uf(1,k) = u(1,k);              %uncomment for not filter the boundaries of the elements
% end

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Energy Based Filter, second try 22-02-2021
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % Differences in energy
% %
% %Legentre
% 
% %By element and mode
% RHSm = Etm-EEnfm-EEm-EEdfm;
% dle = sum(RHSm,1) - Eto;
% 
% % %Lobatto
% % dlo = sum(Etml-EEnfml-EEml-EEdfml,1) - Etol;
% 
% % % %%%
% % % % Lobatto (fails)
% % % %%%
% % % %Correction storage
% % % tol = 1E-8;
% % % corr = zeros(N+1,K);
% % % for k=1:K
% % %    if(abs(dlo(k))>= tol)
% % %        for i=N+1:-1:1
% % % 
% % %            if sign(Etml(i,k))==sign(dlo(k))
% % %                corr(i,k) = corr(i,k) + Etml(i,k);
% % %                if abs(sum(corr(:,k))) >= abs(dlo(k))
% % %                    break
% % %                end
% % %            end   
% % %        end
% % %        corr(i,k) = dlo(k) - sum(corr(i+1:end,k));
% % %    end
% % % end
% % % 
% % % % Correction applied
% % % coa = Etml - EEnfml - corr;
% % % 
% % % % Etml corrected
% % % Etmlc = Etml - corr;
% % %  
% % % % Going back to nodal
% % % % uml = invL*(u.^2);
% % % % Etml = zeros(N+1,K);
% % % % for i = 1:K
% % % %    lode = L.*uml(:,i)';
% % % %    inte = sum(w'.*(lode),1);
% % % %    Etml(:,i) = inte'.*J(:,i);
% % % % end
% % % 
% % % % Recovering modal values
% % % Etmlc2 = Etmlc.*(1./J);
% % % int1 = sum(w'.*L);
% % % umlc = zeros(N+1,K);
% % % for i = 1:K
% % %     umlc(:,i) = Etmlc2(:,i).*(1./int1');
% % % end
% % % 
% % % Etml2 = zeros(N+1,K);
% % % for i = 1:K
% % %    lode2 = L.*umlc(:,i)';
% % %    inte2 = sum(w'.*(lode2),1);
% % %    Etml2(:,i) = inte2'.*J(:,i);
% % % end
% % % 
% % % % Recovering solution filtered
% % % u2c = L*umlc;
% % % u = real(sqrt(u2c)).*sign(u);
% 
% %%%
% % Legendre
% %%%
% %Correction storage
% %Correction storage
% alp = -log(eps);
% umc = um;
% tol = 1E-5;
% corr = zeros(N+1,K);
% 
% %eras could be RHSm, or EEnfm (from what to take to correct)
% % eras = RHSm;
% % eras = -EEnfm;
% % eras = -EEnfm - EEm;
% % eras = EEnfm + EEm;
% % eras = Etm - EEnfm;
% eras = EEnfm+EEm+EEdfm;
% 
% %Energy balance slope in the wholde domain
% Esiw = E(tstep+1)+Ev(tstep+1)+Evf(tstep+1)+Enlf(tstep+1);
% Espw = E(tstep)+Ev(tstep)+Evf(tstep)+Enlf(tstep);
% Ediffw = (Esiw-Espw);
% Eslopew = Ediffw/dt;
% tolw = dt;
% 
% %Energy balance slope by subdomain
% % Esp = E(tstep)+Ev(tstep)+Evf(tstep)+Enlf(tstep);
% % Esi = E(tstep+1)+Ev(tstep+1)+Evf(tstep+1)+Enlf(tstep+1);
% % Esp = E(tstep);
% % Esi = E(tstep+1);
% Esp = EEt(tstep,:)+dEEt(tstep,:)+dfEEt(tstep,:)+nfEEt(tstep,:);
% Esi = EEt(tstep+1,:)+dEEt(tstep+1,:)+dfEEt(tstep+1,:)+nfEEt(tstep+1,:);
% Ediff = (Esi-Esp);
% Eslope = Ediff/dt;
% % tole = 1E-5;
% tole = dt;
% 
% disp([Eslope Eslopew time])
% 
% %changing dle
% % dle = Ediff;
% 
% % if Eslopew > tolw
% %     for k=1:K
% %     %    if(abs(dle(k))>= tol)
% %     %    if(dle(k)>= tol)
% %     %    if(dle(k)> tol)
% %         if (Eslope(k) > tole)
% %            for i=N+1:-1:1
% %                % Computing corr from EEnfm
% %     %            if sign(EEnfm(i,k))==sign(dle(k))
% %     %                corr(i,k) = corr(i,k) + EEnfm(i,k);
% %     %                if abs(sum(corr(:,k))) >= abs(dle(k))
% %     % %                    umc(i:N+1,k) = 0.0;
% %     %                    break
% %     %                end
% %     %            end   
% %                % Computing corr from eras
% %                if sign(eras(i,k))==sign(dle(k))
% %                    corr(i,k) = corr(i,k) + eras(i,k);
% %                    if abs(sum(corr(:,k))) >= abs(dle(k))
% %     %                    umc(i:N+1,k) = 0.0;
% %                        break
% %                    end
% %                end  
% %            end
% %            filterdiag = ones(N+1,1);
% %     %        s = ((N+1)-i)/2;
% % %            s = ((N+1)-i);
% %            s = ((N+1)-i)*100;
% %     %        s = i;
% %     %        Nc = i;
% %     %        Nc = N-4;
% %            for j=Nc:N
% %     % Exponential filter
% %                if Nc==N
% %                    filterdiag(j+1) = 1.0;
% %                else
% %                    filterdiag(j+1) = exp(-alp*((j-Nc)/(N-Nc))^s);
% %                end
% %     % Zeros all
% %     %            if Nc==N
% %     %                filterdiag(j+1) = 0.0;
% %     %            else
% %     %                filterdiag(j+1) = 0.0;
% %     %            end
% %            end
% %     %        i
% %     %        k
% %     %        filterdiag
% % 
% %     %
% %     %  This part is valid only with zeros all option
% %     %
% %     %        corr(i,k) = dle(k) - sum(corr(i+1:end,k));
% %     %        if abs(eras(i,k)) == 0
% %     %            filterdiag(i) = 1.0;
% %     %        else
% %     %            if abs(eras(i,k))<tol
% %     %                filterdiag(i) = 1.0;
% %     %            else
% %     %                filterdiag(i) = 1.0-((abs(eras(i,k)) - abs(corr(i,k)))/abs(eras(i,k)));
% %     % %            filterdiag(i) = ((abs(eras(i,k)) - abs(corr(i,k)))/abs(eras(i,k)));
% %     %            end
% %     %        end
% %     %        %always ensuring the zero polynomial
% %     % %        filterdiag(1) = 1.0;
% %     % %        filterdiag(2) = 1.0;
% %     %        
% %     % %        %always ensuring compatibility of first order polynomial
% %     % %        if filterdiag(2) == 0.0
% %     % %             filterdiag(2) = sqrt(2/3);
% %     % %        end
% % 
% %     %
% %     % End part only valid with zeros all zeros all
% %     %
% % 
% %            %computing new modal values
% %            umc(:,k) = filterdiag.*umc(:,k);
% % 
% %     %        %always ensuring compatibility of first order polynomial
% %     %        if filterdiag(2) == 0.0
% %     %             umc(2,k) = (u(N+1,k)-u(1,k))/(2*V(N+1,2));
% %     %        end
% %        end
% %     end
% %     % pause
% %     % Correction applied
% %     coa = Etm - EEnfm - EEm - EEdfm - corr;
% % 
% %     %Legentre
% %     dle2 = sum(coa,1) - Eto;
% % 
% %     % % Etml corrected
% %     % Etmc = Etm - corr;
% %     %  
% %     % % Going back to nodal
% %     % % Etm =(invV*u).^2.*J(1,:);
% %     % invVu = sqrt(Etmc.*(1./J)).*sign(invV*u);
% % 
% %     % umc(N-4:end,:) = 0.0;
% %     % umc(N:end,:) = 0.5*umc(N:end,:);
% %     % umc(corr~=0) = 0;
% %     % umc(corr~=0) = (1/(1+a))*umc(corr~=0);
% %     % umc = umc + (sqrt(abs(corr)).*sign(corr))*(0.1/(1+a));
% %     % umc = umc + (sqrt(abs(corr)).*sign(corr));
% % %     u = V*umc;
% % 
% %     % %Ensuring continuity
% %     % for k = 2:K
% %     %     u(Np,k-1) = (u(2,k) + u(Np-1,k-1))/2;
% %     %     u(1,k) = u(Np,k-1);
% %     % end
% %     % for k = 2:K
% %     %     u(Np,k-1) = (u(1,k) + u(Np,k-1))/2;
% %     %     u(1,k) = u(Np,k-1);
% %     % end
% % end
% 
% 
% % %Ensuring continuity
% % for k = 2:K
% %     u(Np,k-1) = (u(2,k) + u(Np-1,k-1))/2;
% %     u(1,k) = u(Np,k-1);
% % end
% for k = 2:K
%     u(Np,k-1) = (u(1,k) + u(Np,k-1))/2;
%     u(1,k) = u(Np,k-1);
% end




