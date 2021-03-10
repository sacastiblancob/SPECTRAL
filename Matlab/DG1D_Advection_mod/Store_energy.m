%This script stores the energy calculations (Calc_energy) in the energy
%vectors (Init_energy) for posprocessing comparisons

Globals1D

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% computing energy before filtering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Et = u.^2;

%computing integral with default nodes
EEt(tstep+1,:) = w*(Et.*J);

E(tstep+1) = sum(EEt(tstep+1,:));

%
% advection flux energy total and by mode
%
% Legendre basis
nfEEt(tstep+1,:) = -EEnf;
nfEEtm(:,tstep+1) = reshape(-EEnfm,1,(N+1)*K);

% Lobatto basis
nfEEtml(:,tstep+1) = reshape(-EEnfml,1,(N+1)*K);

%
% Total energy in the system, advection flux
%
%Legendre basis
Enlf(tstep+1) = -1.0*(sum(EEnf));

%    Evtot = sum(EE);
%    Ev(tstep+1) = Ev(tstep) + Evtot;

%Storing times
T(tstep+1) = time;

%Energy per mode
%Legendre basis
um = invV*u;
Etm =(invV*u).^2.*J(1,:);
EEtm(:,tstep+1) = reshape(Etm,1,(N+1)*K);

%Lobato basis
% Etml =(invL*u).^2.*J(1,:); %This does not work
uml1 = invL*u;
uml = invL*(u.^2);
Etml = zeros(N+1,K);
for i = 1:K
   lode = L.*uml(:,i)';
   inte = sum(w'.*(lode),1);
   Etml(:,i) = inte'.*J(:,i);
end
EEtml(:,tstep+1) = reshape(Etml,1,(N+1)*K);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Differences in energy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Legentre
% A = sum(Etm-EEnfm,1)
% B = sum(Etm,1)
% sum(EEnfm,1)
RHSm = Etm-EEnfm;
dle = sum(RHSm,1) - Eto

%Lobatto
dlo = sum(Etml-EEnfml,1) - Etol;

% % %%%
% % % Lobatto (fails)
% % %%%
% % %Correction storage
% % tol = 1E-8;
% % corr = zeros(N+1,K);
% % for k=1:K
% %    if(abs(dlo(k))>= tol)
% %        for i=N+1:-1:1
% % 
% %            if sign(Etml(i,k))==sign(dlo(k))
% %                corr(i,k) = corr(i,k) + Etml(i,k);
% %                if abs(sum(corr(:,k))) >= abs(dlo(k))
% %                    break
% %                end
% %            end   
% %        end
% %        corr(i,k) = dlo(k) - sum(corr(i+1:end,k));
% %    end
% % end
% % 
% % % Correction applied
% % coa = Etml - EEnfml - corr;
% % 
% % % Etml corrected
% % Etmlc = Etml - corr;
% %  
% % % Going back to nodal
% % % uml = invL*(u.^2);
% % % Etml = zeros(N+1,K);
% % % for i = 1:K
% % %    lode = L.*uml(:,i)';
% % %    inte = sum(w'.*(lode),1);
% % %    Etml(:,i) = inte'.*J(:,i);
% % % end
% % 
% % % Recovering modal values
% % Etmlc2 = Etmlc.*(1./J);
% % int1 = sum(w'.*L);
% % umlc = zeros(N+1,K);
% % for i = 1:K
% %     umlc(:,i) = Etmlc2(:,i).*(1./int1');
% % end
% % 
% % Etml2 = zeros(N+1,K);
% % for i = 1:K
% %    lode2 = L.*umlc(:,i)';
% %    inte2 = sum(w'.*(lode2),1);
% %    Etml2(:,i) = inte2'.*J(:,i);
% % end
% % 
% % % Recovering solution filtered
% % u2c = L*umlc;
% % u = real(sqrt(u2c)).*sign(u);

%%%
% Legendre
%%%
%Correction storage
alp = -log(eps);
umc = um;
tol = 1E-8;
corr = zeros(N+1,K);
for k=1:K
%    if(abs(dle(k))>= tol)
%    if(dle(k)>= tol)
   if(dle(k)> 0)
       for i=N+1:-1:1
           % Computing corr from EEnfm
%            if sign(EEnfm(i,k))==sign(dle(k))
%                corr(i,k) = corr(i,k) + EEnfm(i,k);
%                if abs(sum(corr(:,k))) >= abs(dle(k))
% %                    umc(i:N+1,k) = 0.0;
%                    break
%                end
%            end
           % Computing corr directly from the sum
           if sign(RHSm(i,k))==sign(dle(k))
               corr(i,k) = corr(i,k) + EEnfm(i,k);
               if abs(sum(corr(:,k))) >= abs(dle(k))
%                    umc(i:N+1,k) = 0.0;
                   break
               end
           end  
       end
       filterdiag = ones(N+1,1);
       s = ((N+1)-i)/2;
%        s = ((N+1)-i);
%        s = i;
       Nc = i;
       for j=Nc:N
% Exponential filter
           if Nc==N
               filterdiag(j+1) = 1.0;
           else
               filterdiag(j+1) = exp(-alp*((j-Nc)/(N-Nc))^s);
           end
% Zeros all
%            if Nc==N
%                filterdiag(j+1) = 1.0;
%            else
%                filterdiag(j+1) = 0.0;
%            end
       end
%        i
%        k
%        filterdiag
       umc(:,k) = filterdiag.*umc(:,k);
       corr(i,k) = dle(k) - sum(corr(i+1:end,k));
   end
end
% pause
% Correction applied
coa = Etm - EEnfm - corr;

%Legentre
dle2 = sum(coa,1) - Eto

% % Etml corrected
% Etmc = Etm - corr;
%  
% % Going back to nodal
% % Etm =(invV*u).^2.*J(1,:);
% invVu = sqrt(Etmc.*(1./J)).*sign(invV*u);

% umc(N-4:end,:) = 0.0;
% umc(N:end,:) = 0.5*umc(N:end,:);
% umc(corr~=0) = 0;
% umc(corr~=0) = (1/(1+a))*umc(corr~=0);
% umc = umc + (sqrt(abs(corr)).*sign(corr))*(0.1/(1+a));
% umc = umc + (sqrt(abs(corr)).*sign(corr));
u = V*umc;



% %Ensuring continuity
% for k = 2:K
%     u(Np,k-1) = (u(2,k) + u(Np-1,k-1))/2;
%     u(1,k) = u(Np,k-1);
% end
% for k = 2:K
%     u(Np,k-1) = (u(1,k) + u(Np,k-1))/2;
%     u(1,k) = u(Np,k-1);
% end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% storing energy for Hovmoller plots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% First term of energy equation u^2 (actual energy of the solution)
EH(tstep,:) = reshape(u.^2,1,K*(N+1));

% Second term of energy equation -a*u^2 (advection flux)
EHnlf(tstep,:) = reshape((2/3)*u.^3,1,K*(N+1));

% Coordinates
XCoor(tstep,:) = reshape(x,1,K*(N+1));

% Times
Times(tstep,:) = time*ones(1,K*(N+1));

%%%
%end
%%%


