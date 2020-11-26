%This script stores the energy calculations (Calc_energy) in the energy
%vectors (Init_energy) for posprocessing comparisons

Globals1D

% computing energy before filtering
if iint == 0
    Et = u.^2;

    EEt(tstep+1,:) = w*(Et.*J);    %computing integral with default nodes
else
    umd = [invV*u ; zeros(length(rd)-length(r),Elements)];    %u modal and zero-padding
    ud = Vd*umd;    %to nodal with more modes
    Et = ud.^2;     %squaring

    EEt(tstep+1,:) = wd*(Et.*Jd);   %computing integral with more nodes
    umd = umd*0.0;
    ud = ud*0.0;
end

E(tstep+1) = sum(EEt(tstep+1,:));

%     %Computing dissipate energy due to viscosity
%     %duv = (1./J).*Dr*u;
%     dus = (Dr*u./J).^2;
%     %function of N, Elements, CFL, xL and xR
%     %Evele = 2*0.008059650327027*epsilon*w*(dus.*J);
%     Evele = dt*epsilon*w*(dus.*J);
%     %Evele = epsilon*w*(dus.*J);
dEEt(tstep+1,:) = -EE;
dEEtm(:,:,tstep+1) = -EEm;
dfEEt(tstep+1,:) = -EEdf;
dfEEtm(:,:,tstep+1) = - EEdfm;
nfEEt(tstep+1,:) = -EEnf;
nfEEtm(:,:,tstep+1) = - EEnfm;
Ev(tstep+1) = -1.0*(sum(EE));
Evf(tstep+1) = -1.0*(sum(EEdf));
Enlf(tstep+1) = -1.0*(sum(EEnf));

%    Evtot = sum(EE);
%    Ev(tstep+1) = Ev(tstep) + Evtot;

%Storing times
T(tstep+1) = time;

%Energy per mode
Etm =(invV*u).^2.*J(1,:);
EEtm(:,:,tstep+1) = Etm;

% %
% %Energy Based Filter
% %
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







