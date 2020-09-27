

Globals1D

%Plot of global energy behaviour

% hold on
% plot(T,E)
% hold on
% plot(T,Ev)
% hold on
% plot(T,Evf)
% hold on
% plot(T,Enlf)
% hold on
% Ebal = E+Ev+Evf+Enlf;
% plot(T,Ebal)
% legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
% %legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','E1','E_{vd}1','E_{vf}1','E_{nlf}1','E+E_{vd}+E_{vf}+E_{nlf}1','Location','eastoutside')
% title('Energy vs Time, \nu=0.0, \Omega = whole domain')
% %ylim([-E(1)/10 E(1)+E(1)/5])
% ylim([-0.1 1.6])
% xlim([T(1) T(tstep+1)])
% %xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')
% 
% a=2;
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
% title('Energy vs Time, \nu=0.0, \Omega = 5th domain')
% %ylim([liminf limsup])
% %ylim([-0.22 0.32])
% ylim([-0.5 0.2])
% xlim([T(1) T(tstep+1)])
% %xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')







