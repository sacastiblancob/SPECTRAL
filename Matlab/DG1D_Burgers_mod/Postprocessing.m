

Globals1D

%Plot of global energy behaviour

hold on
plot(T,E)
hold on
plot(T,Ev)
hold on
plot(T,Evf)
hold on
plot(T,Enlf)
hold on
Ebal = E+Ev+Evf+Enlf;
Eotro = E+Ev;
plot(T,Eotro)
legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
%legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','E1','E_{vd}1','E_{vf}1','E_{nlf}1','E+E_{vd}+E_{vf}+E_{nlf}1','Location','eastoutside')
title('Energy vs Time, \nu=0.0, \Omega = whole domain')
%ylim([-E(1)/10 E(1)+E(1)/5])
ylim([-0.1 1.4])
%xlim([T(1) T(tstep+1)])
xlim([T(1) 0.8])
ylabel('Energy')
xlabel('Time')
% 

% a=6;
% liminf = min([EEt(:,a);dEEt(:,a);dfEEt(:,a);nfEEt(:,a)])-EEt(1,a)/10;
% limsup = max([EEt(:,a);dEEt(:,a);dfEEt(:,a);nfEEt(:,a)])+EEt(1,a)/1.5;
% 
% plot(T,EEt(:,a))
% hold on
% plot(T,dEEt(:,a))
% plot(T,dfEEt(:,a))
% plot(T,nfEEt(:,a))
% plot(T,Econ(:,a))
% legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
% title(strcat('Energy vs Time, \nu=0.0, \Omega = ',string(a),'th domain'))
% %ylim([liminf limsup])
% %ylim([-0.22 0.32])
% ylim([-0.5 0.5])
% xlim([T(1) T(tstep+1)])
% %xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')

% %%% Proofs
% eo = EEt(1,:);
% eom = EEtm(:,:,1);
% pt = 1000;
% e_sum = EEt(pt,:) + dEEt(pt,:) + dfEEt(pt,:) + nfEEt(pt,:);
% e_summ = EEtm(:,:,pt) + dEEtm(:,:,pt) + dfEEtm(:,:,pt) + nfEEtm(:,:,pt);
% e_m = sum(e_summ);








