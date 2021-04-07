

Globals1D

% Econ = EEt + dEEt + dfEEt + nfEEt;
% %Econ = EEt + dEEt + nfEEt;
% 
% % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % % Plot of global energy behaviour
% % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% figure(2)
% 
% hold on
% plot(T,E,'LineWidth',1.5)
% hold on
% plot(T,Ev,'LineWidth',1.5)
% hold on
% plot(T,Evf,'LineWidth',1.5)
% hold on
% plot(T,Enlf,'LineWidth',1.5)
% hold on
% Ebal = E+Ev+Evf+Enlf;
% % Eotro = E+Ev;
% Eotro = Ebal;
% plot(T,Eotro,'LineWidth',2)
% legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
% legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','E1','E_{vd}1','E_{vf}1','E_{nlf}1','E+E_{vd}+E_{vf}+E_{nlf}1','Location','eastoutside')
% title(strcat('Energy vs Time, \nu=',string(epsilon),', \Omega = whole domain'))
% ylim([-E(1)/10 E(1)+E(1)/5])
% % ylim([-0.9 1.4])
% xlim([T(1) T(tstep+1)])
% % xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')
% grid on
% 
% %Slope of Ebal
% EbalEsl = zeros(Nsteps+1,1);
% for i=1:Nsteps
%     EbalEsl = (Ebal(i+1)-Ebal(i))/dt;
% end
% 
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % Plot of global energy behaviour by mode a
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% figure(3)
% 
% ap = 1;
% a = ap;
% 
% liminf = min([EEt(:,a);dEEt(:,a);dfEEt(:,a);nfEEt(:,a)])-EEt(1,a)/10;
% limsup = max([EEt(:,a);dEEt(:,a);dfEEt(:,a);nfEEt(:,a)])+EEt(1,a)/1.5;
% 
% plot(T,EEt(:,a),'LineWidth',1.5)
% hold on
% plot(T,dEEt(:,a),'LineWidth',1.5)
% plot(T,dfEEt(:,a),'LineWidth',1.5)
% plot(T,nfEEt(:,a),'LineWidth',1.5)
% plot(T,Econ(:,a),'LineWidth',1.5)
% legend('E','E_{vd}','E_{vf}','E_{nlf}','E+E_{vd}+E_{vf}+E_{nlf}','Location','eastoutside')
% title(strcat('Energy vs Time, \nu =',string(epsilon),', \Omega = ',string(a),'th domain'))
% %ylim([liminf limsup])
% % ylim([-0.22 0.32])
% ylim([-0.3 0.4])
% % ylim([min([EEt(:,a) ; dEEt(:,a); dfEEt(:,a); nfEEt(:,a); Econ(:,a)])-0.05, max([EEt(:,a) ; dEEt(:,a); dfEEt(:,a); nfEEt(:,a); Econ(:,a)])+0.05])
% xlim([T(1) T(tstep+1)])
% %xlim([T(1) 0.8])
% ylabel('Energy')
% xlabel('Time')
% grid on
% 
% %%% Proofs
% eo = EEt(1,:);
% eom = EEtm(:,:,1);
% pt = 1000;
% e_sum = EEt(pt,:) + dEEt(pt,:) + dfEEt(pt,:) + nfEEt(pt,:);
% e_summ = EEtm(:,pt) + dEEtm(:,pt) + dfEEtm(:,pt) + nfEEtm(:,pt);
% e_m = sum(e_summ);
% 
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %Hovmoller plots
% % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % 
% % %
% % % ENERGY TERMS
% % %
% % HS = EH + EHv;
% % 
% % figure(4)
% % 
% % subplot(2,3,1)
% % surf(XCoor,Times,EH)
% % xlabel('x')
% % ylabel('t')
% % title('u^{2}')
% % axis([xL xR 0 FinalTime]);
% % view(0,90);
% % shading interp
% % colorbar
% % 
% % subplot(2,3,2)
% % surf(XCoor,Times,EHnlf)
% % xlabel('x')
% % ylabel('t')
% % title('-(2/3)u^3')
% % axis([xL xR 0 FinalTime]);
% % view(0,90);
% % shading interp
% % colorbar
% % 
% % subplot(2,3,3)
% % surf(XCoor,Times,EHvf)
% % xlabel('x')
% % ylabel('t')
% % title('2 \nu (du/dx)')
% % % axis([xL xR 0 FinalTime -0.05 0.05]);
% % view(0,90);
% % shading interp
% % colorbar
% % % caxis([-0.05 0.05])
% % 
% % subplot(2,3,4)
% % surf(XCoor,Times,EHv)
% % xlabel('x')
% % ylabel('t')
% % title('(du/dx)^2')
% % % axis([xL xR 0 FinalTime 0 1]);
% % view(0,90);
% % shading interp
% % colorbar
% % % caxis([0 5])
% % 
% % subplot(2,3,5)
% % surf(XCoor,Times,HS)
% % xlabel('x')
% % ylabel('t')
% % title('u^2 - 2 \nu (du/dx)^2')
% % % axis([xL xR 0 FinalTime 0 1]);
% % view(0,90);
% % shading interp
% % colorbar
% % % caxis([0 5])
% 
% XCoor = kron((xL+(xR-xL)/(2*(N+1)*K):(xR-xL)/((N+1)*K):xR)',ones(1,tstep));
% 
% %
% % MODAL ENERGY
% %
% EMS = EEtm + dEEtm + dfEEtm + nfEEtm;
% 
% %Energy balance by mode
% EBM = zeros(K,tstep+1);
% for i=1:K
%     EBM(i,:) = sum(EMS((i-1)*(N+1)+1:i*(N+1),:),1);
% end
% 
% dtpl = floor(tstep/40);
% XPL = XCoor;
% XPL = XPL(:,1:dtpl:tstep);
% YPL = Times';
% YPL = YPL(:,1:dtpl:tstep);
% 
% % % figure(5)
% % % 
% % % subplot(2,3,1)
% % % ZPL = EEtm(:,1:tstep);
% % % ZPL = ZPL(:,1:dtpl:tstep);
% % % surf(XPL,YPL,ZPL)
% % % xlabel('x')
% % % ylabel('t')
% % % title('Actual energy of the system')
% % % % axis([xL xR 0 FinalTime]);
% % % view(0,90);
% % % xticks(xL:(xR-xL)/K:xR)
% % % % shading interp
% % % colorbar
% % % 
% % % subplot(2,3,2)
% % % ZPL = dEEtm(:,1:tstep);
% % % ZPL = ZPL(:,1:dtpl:tstep);
% % % surf(XPL,YPL,ZPL)
% % % xlabel('x')
% % % ylabel('t')
% % % title('Viscous energy')
% % % % axis([xL xR 0 FinalTime]);
% % % view(0,90);
% % % xticks(xL:(xR-xL)/K:xR)
% % % % shading interp
% % % colorbar
% % % 
% % % subplot(2,3,3)
% % % ZPL = dfEEtm(:,1:tstep);
% % % ZPL = ZPL(:,1:dtpl:tstep);
% % % surf(XPL,YPL,ZPL)
% % % xlabel('x')
% % % ylabel('t')
% % % title('Viscous fluxes')
% % % % axis([xL xR 0 FinalTime -0.05 0.05]);
% % % view(0,90);
% % % xticks(xL:(xR-xL)/K:xR)
% % % % shading interp
% % % colorbar
% % % % caxis([-0.05 0.05])
% % % 
% % % subplot(2,3,4)
% % % ZPL = nfEEtm(:,1:tstep);
% % % ZPL = ZPL(:,1:dtpl:tstep);
% % % surf(XPL,YPL,ZPL)
% % % xlabel('x')
% % % ylabel('t')
% % % title('Non-linear fluxes')
% % % % axis([xL xR 0 FinalTime 0 1]);
% % % view(0,90);
% % % xticks(xL:(xR-xL)/K:xR)
% % % % shading interp
% % % colorbar
% % % % caxis([0 5])
% % % 
% % % subplot(2,3,5)
% % % ZPL = EMS(:,1:tstep);
% % % ZPL = ZPL(:,1:dtpl:tstep);
% % % surf(XPL,YPL,ZPL)
% % % xlabel('x')
% % % ylabel('t')
% % % title('Sum of energy per Mode')
% % % % axis([xL xR 0 FinalTime 0 1]);
% % % view(0,90);
% % % xticks(xL:(xR-xL)/K:xR)
% % % % shading interp
% % % colorbar
% % % % caxis([0 5])
% % % 
% % % % sumEMS = [sum(EMS,1);sum(EMS,1)];
% % % subplot(2,3,6)
% % % XPL = kron((1:K)',ones(1,tstep));
% % % XPL = XPL(:,1:dtpl:tstep);
% % % YPL = Times(:,1:K)';
% % % YPL = YPL(:,1:dtpl:tstep);
% % % ZPL = EBM(:,1:tstep);
% % % ZPL = ZPL(:,1:dtpl:tstep);
% % % surf(XPL,YPL,ZPL)
% % % xlabel('x')
% % % ylabel('t')
% % % title('Sum of energy per Element')
% % % % axis([xL xR 0 FinalTime 0 1]);
% % % view(0,90);
% % % % shading interp
% % % colorbar
% % % % caxis([0 5])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Legendre basis
% figure(7)
% 
% % for m = 1:N+1
% for m = 1:N
%     %Plotting the actual energy by mode for the element ap
%     subplot(2,2,1)
%     plot3(T,ones(length(T),1)*(m-1),EEtm(((ap-1)*(N+1))+m,:),'Color','b')
%     grid on
%     hold on
%     
%     %Plotting the energy advection fluxes by mode for the element ap
%     subplot(2,2,4)
%     if m==12
%         plot3(T,ones(length(T),1)*(m-1),nfEEtm(((ap-1)*(N+1))+m,:),'Color',[255/255 100/255 0])
%     else
%         plot3(T,ones(length(T),1)*(m-1),nfEEtm(((ap-1)*(N+1))+m,:),'Color','m')
%     end
%     grid on
%     hold on
%     
%     %Plotting the energy advection fluxes by mode for the element ap
%     subplot(2,2,2)
%     if m==12
%         plot3(T,ones(length(T),1)*(m-1),dEEtm(((ap-1)*(N+1))+m,:),'Color','b')
%     else
%         plot3(T,ones(length(T),1)*(m-1),dEEtm(((ap-1)*(N+1))+m,:),'Color','r')
%     end
%     grid on
%     hold on
%     
%     %Plotting the energy advection fluxes by mode for the element ap
%     subplot(2,2,3)
%     plot3(T,ones(length(T),1)*(m-1),dfEEtm(((ap-1)*(N+1))+m,:),'Color',[255/255 150/255 0])
%     grid on
%     hold on
% 
% %     %Plotting the sum of energies by mode for the element ap
% %     subplot(2,2,3)
% %     plot3(T,ones(length(T),1)*(m-1),EMS(((ap-1)*(N+1))+m,:),'Color','g')
% %     grid on
% %     hold on
%     
% end
% 
% %Plotting the actual energy by mode for the element ap
% subplot(2,2,1)
% xlabel('Time')
% ylabel('Mode')
% zlabel('Energy')
% set(gca, 'Ydir', 'reverse')
% sumle = sum(EEtm(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
% plot3(T,ones(length(T),1)*(-1),sumle,'Color','k')
% grid on
% zlim([-0.1 0.3])
% ylim([-1.1 K+0.1])
% xlim([0 time])
% % title(strcat('Actual Energy by mode of ',string(ap),'th domain (Black is sum)'))
% title('E')
% 
% %Plotting the energy advection fluxes by mode for the element ap
% subplot(2,2,4)
% xlabel('Time')
% ylabel('Mode')
% zlabel('Energy')
% set(gca, 'Ydir', 'reverse')
% sumlef = sum(nfEEtm(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
% plot3(T,ones(length(T),1)*(-1),sumlef,'Color','k')
% zlim([-0.3 0.3])
% ylim([-1.1 K+0.1])
% xlim([0 time])
% grid on
% % title(strcat('Energy Advection Flux by mode of ',string(ap),'th domain (Black is sum)'))
% title('E_{nlf}')
% 
% %Plotting the energy advection fluxes by mode for the element ap
% subplot(2,2,2)
% xlabel('Time')
% ylabel('Mode')
% zlabel('Energy')
% set(gca, 'Ydir', 'reverse')
% sumlef = sum(dEEtm(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
% plot3(T,ones(length(T),1)*(-1),sumlef,'Color','k')
% zlim([-0.02 0.08])
% ylim([-1.1 K+0.1])
% xlim([0 time])
% grid on
% % title(strcat('Energy Advection Flux by mode of ',string(ap),'th domain (Black is sum)'))
% title('E_{vd}')
% 
% %Plotting the energy advection fluxes by mode for the element ap
% subplot(2,2,3)
% xlabel('Time')
% ylabel('Mode')
% zlabel('Energy')
% set(gca, 'Ydir', 'reverse')
% sumlef = sum(dfEEtm(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
% plot3(T,ones(length(T),1)*(-1),sumlef,'Color','k')
% zlim([-0.03 0.03])
% ylim([-1.1 K+0.1])
% xlim([0 time])
% grid on
% % title(strcat('Energy Advection Flux by mode of ',string(ap),'th domain (Black is sum)'))
% title('E_{vf}')
% sgtitle(strcat('\Omega =',string(ap),'th subdomain'))
% 
% % % %Plotting the sum of energies by mode for the element ap
% % % subplot(2,2,3)
% % % xlabel('Time')
% % % ylabel('Mode')
% % % zlabel('Energy')
% % % set(gca, 'Ydir', 'reverse')
% % % sumlet = sum(EMS(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
% % % plot3(T,ones(length(T),1)*(-1),sumlet,'Color','k')
% % % grid on
% % % % title(strcat('Sum of Energies by mode of ',string(ap),'th domain (Black is sum)'))
% % % title(strcat('Sum of Energies by mode: ',string(ap),'th domain'))
% 
% % subplot(2,2,4)
% % for i = 1 : Elements
% % %         plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
% % %           'Markersize', 1, 'LineWidth', 2.0, 'Color','b');
% %     plot ( x(:,i), u(:,i), shapestr{1+rem(i,2)}, ...
% %       'Markersize', 1, 'LineWidth', 1.5);
% %     hold all
% % end
% % title('Solution')
% % xlabel('X')
% % ylabel('U')
% % grid on
% 
% % % %%%
% % % % plot of u for desired t
% % % %%%
% % % 
% % % uvec = reshape(u,[1 (N+1)*K]);      %stable
% % % uunsv = reshape(uuns,[1 (N+1)*K]);  %unstable
% % % xvec = reshape(x,[1 (N+1)*K]);
% % % 
% % % plot(xvec,uvec,'LineWidth',1.5)
% % % hold on
% % % plot(xvec,uunsv,'LineWidth',1.5)
% % % legend('\nu=5E^{-2}','\nu=6E^{-5}')
% % % xlabel('X')
% % % ylabel('U')
% % % ylim([-1.2 1.2])
% % % title('Solution at t=0.4s')
% % % grid on

%%%%%%%

Econa = EEta + dEEta + dfEEta + nfEEta;
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % Plot of global energy behaviour ANALYTICAL
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(8)

hold on
plot(T,Ea,'LineWidth',1.5)
hold on
plot(T,Eva,'LineWidth',1.5)
hold on
plot(T,Evfa,'LineWidth',1.5)
hold on
plot(T,Enlfa,'LineWidth',1.5)
hold on
Ebala = Ea+Eva+Evfa+Enlfa;
% Eotro = E+Ev;
Eotroa = Ebala;
plot(T,Eotroa,'LineWidth',2)
legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','Location','eastoutside')
% legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','E1 a','E_{vd}1 a','E_{vf}1 a','E_{nlf}1 a','E+E_{vd}+E_{vf}+E_{nlf}1 a','Location','eastoutside')
title(strcat('Energy vs Time, \nu=',string(epsilon),', \Omega = whole domain (analytical)'))
ylim([-E(1)/10 E(1)+E(1)/5])
% ylim([-0.9 1.4])
xlim([T(1) T(tstep+1)])
% xlim([T(1) 0.8])
ylabel('Energy (analytical)')
xlabel('Time')
grid on

% %Slope of Ebal
% EbalEsl = zeros(Nsteps+1,1);
% for i=1:Nsteps
%     EbalEsl = (Ebal(i+1)-Ebal(i))/dt;
% end

% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Plot of global energy behaviour by mode a ANALYTICAL
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(9)

ap = 4;
a = ap;

liminf = min([EEta(:,a);dEEta(:,a);dfEEta(:,a);nfEEta(:,a)])-EEta(1,a)/10;
limsup = max([EEta(:,a);dEEta(:,a);dfEEta(:,a);nfEEta(:,a)])+EEta(1,a)/1.5;

plot(T,EEta(:,a),'LineWidth',1.5)
hold on
plot(T,dEEta(:,a),'LineWidth',1.5)
plot(T,dfEEta(:,a),'LineWidth',1.5)
plot(T,nfEEta(:,a),'LineWidth',1.5)
plot(T,Econa(:,a),'LineWidth',1.5)
legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','Location','eastoutside')
title(strcat('Energy vs Time, \nu =',string(epsilon),', \Omega = ',string(a),'th domain (analytical)'))
%ylim([liminf limsup])
% ylim([-0.22 0.32])
ylim([-0.3 0.4])
% ylim([min([EEt(:,a) ; dEEt(:,a); dfEEt(:,a); nfEEt(:,a); Econ(:,a)])-0.05, max([EEt(:,a) ; dEEt(:,a); dfEEt(:,a); nfEEt(:,a); Econ(:,a)])+0.05])
xlim([T(1) T(tstep+1)])
%xlim([T(1) 0.8])
ylabel('Energy  (analytical)')
xlabel('Time')
grid on


% By computing Eotroa and Eotroaint, which is the balance over analytical
% with dealeasing 3/2 rule, the norm of the result over all the time is
% about 1E-2 for two norm, and 1E-4 for Inf norm. Which is insignificant.
% plot(T,Eotroa,'LineWidth',2)
% hold on
% plot(T,Eotroaint,'LineWidth',2)

%%%%%%

figure(10)

hold on
plot(T,Ea,'LineWidth',1.5)
hold on
plot(T,Eva,'LineWidth',1.5)
hold on
plot(T,Evfa,'LineWidth',1.5)
hold on
plot(T,Enlfa,'LineWidth',1.5)
hold on
Ebala = Ea+Evfa+Enlfa;
% Eotro = E+Ev;
Eotroa = Ebala;
plot(T,Eotroa,'LineWidth',2)
legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vf}+E_{nlf} a','Location','eastoutside')
% legend('E a','E_{vd} a','E_{vf} a','E_{nlf} a','E+E_{vd}+E_{vf}+E_{nlf} a','E1 a','E_{vd}1 a','E_{vf}1 a','E_{nlf}1 a','E+E_{vd}+E_{vf}+E_{nlf}1 a','Location','eastoutside')
title(strcat('Energy vs Time, \nu=',string(epsilon),', \Omega = whole domain (analytical)'))
ylim([-E(1)/10 E(1)+E(1)/5])
% ylim([-0.9 1.4])
xlim([T(1) T(tstep+1)])
% xlim([T(1) 0.8])
ylabel('Energy (analytical)')
xlabel('Time')
grid on




