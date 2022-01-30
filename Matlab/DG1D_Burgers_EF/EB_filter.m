%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sergio A. Castiblanco B.
% Pontificia Universidad Javeriana - Bogota
% Hydrosystems Master Program
%
% Energy Balance Based Filter Methodology
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Energy Based Filter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Differences in energy
%
%Legentre

Globals1D

%By element and mode
RHSm = Etm-EEnfm-EEm-EEdfm;
dle = sum(RHSm,1) - Eto;

%Correction storage
alp = -log(eps);
umc = um;
tol = 1E-5;
corr = zeros(N+1,K);

%eras could be RHSm, or EEnfm (from what to take to correct)
% eras = RHSm;
% eras = -EEnfm;
% eras = -EEnfm - EEm;
% eras = EEnfm + EEm;
% eras = Etm - EEnfm;
% eras = EEnfm+EEm+EEdfm;
% eras = -Etm + EEm + EEdfm;
eras = -RHSm;
% eras = EEm + EEdfm;

%Energy balance slope in the wholde domain
Esiw = E(tstep+1)+Ev(tstep+1)+Evf(tstep+1)+Enlf(tstep+1);
Espw = E(tstep)+Ev(tstep)+Evf(tstep)+Enlf(tstep);
Ediffw = (Esiw-Espw);
Eslopew = Ediffw/dt;
tolw = dt;
% tolw = 0;

%Energy balance slope by subdomain
% Esp = E(tstep)+Ev(tstep)+Evf(tstep)+Enlf(tstep);
% Esi = E(tstep+1)+Ev(tstep+1)+Evf(tstep+1)+Enlf(tstep+1);
% Esp = E(tstep);
% Esi = E(tstep+1);
Esp = EEt(tstep,:)+dEEt(tstep,:)+dfEEt(tstep,:)+nfEEt(tstep,:);
Esi = EEt(tstep+1,:)+dEEt(tstep+1,:)+dfEEt(tstep+1,:)+nfEEt(tstep+1,:);
Ediff = (Esi-Esp);
Eslope = Ediff/dt;
% tole = dt/10;
tole = dt;

% disp([Eslope Eslopew time])

%changing dle
dle = Ediff;

if Eslopew > Eslopew0
    for k=1:K
      if (Eslope(k) > Eslope0(k))
            knaf(k) = 0;
           %Searching modes for recovering Energy
           for i=N+1:-1:1
               % Computing corr from eras
               if sign(eras(i,k))==sign(dle(k))
                   corr(i,k) = corr(i,k) + eras(i,k);
                   if abs(sum(corr(:,k))) >= abs(dle(k))
    %                    umc(i:N+1,k) = 0.0;
                       break
                   end
               end  
           end
           filterdiag = ones(N+1,1);
           i = i-1;
           Nc = i;
           if Nc<Ncs(k)
             Ncs(k) = Nc;                     %Computed cutoff
             if filter_ebf == 1
    %            s = ((N+1)-i)/(2*pi);    %Order for Exponential Filter
    %            s = round((Nc)/2)*2;
    %            s = round(Nc/2);
    %            s = ((N+1)-i);
    %            s = ((N+1)-i)*100;
    %            s = i;
%                  s = CFL*((N+1)-i);
                 s = round((N-Nc)/2)*2;
             else
                 s = round((N-Nc)/2)*2;             %Order for Vandeven kind filter (%naaays)
    %            s = N;
             end

           for j=Nc:N
               if Nc==N
    %                    filterdiag(j+1) = 1.0;
                   Fdiags(j+1,k) = 0.0;
                   Fdiags(j,k) = 0.0;
               else
                   ex = (j-Nc)/(N-Nc);
                   if filter_ebf == 1
% %                    Exponential
%                        Fdiags(j+1,k) = exp(-alp*(ex)^s);
                       Fdiags(j+1,:) = min(Fdiags(j+1,:),ones(1,K)*exp(-alp*(ex)^s));
                   else

%     %                    %Vandeven
                       yv = -((ex^(-s))*((-(ex - 1)*ex)^s)*hypergeom([1-s s],s+1,1-ex));
                       yv = yv/s;
%                        Fdiags(j+1,k) = 0 - (factorial(2*s - 1)/(factorial(s - 1)^2))*yv;
                       Fdiags(j+1,:) = min(Fdiags(j+1,:),ones(1,K)*(0 - (factorial(2*s - 1)/(factorial(s - 1)^2))*yv));
                   end
               end
           end
           end

           if filter_ebf == 2
               % Removing Nan's if Vandeven
               Fdiags(isnan(Fdiags)) = 1.0;
           end

       end
    end
        %Exponential
%     for j0=Nc0:N
%        Fdiags(j0+1,logical(knaf)) = ones(1,sum(knaf))*exp(-alp*((j0-Nc0)/(N-Nc0))^(s0));
%     end

% %     %Vandeven
%     for j0=Nc0+1:N
%        ex0 = (j0-Nc0)/(N-Nc0);
%        yv = -((ex0^(-s0))*((-(ex0 - 1)*ex0)^s0)*hypergeom([1-s0 s0],s0+1,1-ex0));
%        yv = yv/s0;
%        Fdiags(j0+1,logical(knaf)) = ones(1,sum(knaf))*(0 - (factorial(2*s0 - 1)/(factorial(s0 - 1)^2))*yv);
%     end
end
umc = Fdiags.*umc;
u = V*umc;

Eslopew0 = CFL*Eslopew;
Eslope0 = CFL*Eslope;

