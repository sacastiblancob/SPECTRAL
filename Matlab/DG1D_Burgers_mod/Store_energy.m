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
Etom =(invV*u).^2.*J(1,:);
EEtm(:,:,tstep+1) = Etom;








