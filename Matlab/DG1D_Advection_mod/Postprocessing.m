

Globals1D

% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % Plot of global energy behaviour
% % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(2)

hold on
plot(T,E)
hold on
plot(T,Enlf)
hold on
Ebal = E+Enlf;
plot(T,Ebal)
legend('E','E_{adv}','E+E_{adv}','Location','eastoutside')
title('Energy vs Time, \nu=0.0, \Omega = whole domain')
ylim([-E(1)/10 E(1)+E(1)/5])
ylim([-0.1 max([E;Enlf;Ebal])+0.1])
xlim([T(1) T(tstep+1)])
xlim([T(1) 0.8])
ylabel('Energy')
xlabel('Time')
grid on

% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Plot of global energy behaviour by mode a
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure(3)

ap=4;
liminf = min([EEt(:,ap);nfEEt(:,ap)])-EEt(1,ap)/10;
limsup = max([EEt(:,ap);nfEEt(:,ap)])+EEt(1,ap)/1.5;

plot(T,EEt(:,ap))
hold on
plot(T,nfEEt(:,ap))
plot(T,Econ(:,ap))
legend('E','E_{adv}','E+E_{adv}','Location','eastoutside')
title(strcat('Energy vs Time, \Omega = ',string(ap),'th domain'))
%ylim([liminf limsup])
%ylim([-0.22 0.32])
%ylim([-0.5 0.5])
ylim([min([EEt(:,ap) ; nfEEt(:,ap); Econ(:,ap)])-0.1, max([EEt(:,ap) ; nfEEt(:,ap); Econ(:,ap)])+0.1])
xlim([T(1) T(tstep+1)])
%xlim([T(1) 0.8])
ylabel('Energy')
xlabel('Time')

grid on

% % %%% Proofs
% % eo = EEt(1,:);
% % eom = EEtm(:,:,1);
% % pt = 1000;
% % e_sum = EEt(pt,:) + dEEt(pt,:) + dfEEt(pt,:) + nfEEt(pt,:);
% % e_summ = EEtm(:,pt) + dEEtm(:,pt) + dfEEtm(:,pt) + nfEEtm(:,pt);
% % e_m = sum(e_summ);
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %Hovmoller plots
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
XCoor = kron((xL+(xR-xL)/(2*(N+1)*K):(xR-xL)/((N+1)*K):xR)',ones(1,tstep));

%
% MODAL ENERGY
%
%Legendre basis
EMS = EEtm + nfEEtm;

%Energy balance by mode
EBM = zeros(K,tstep+1);
for i=1:K
    EBM(i,:) = sum(EMS((i-1)*(N+1)+1:i*(N+1),:),1);
end

%Lobatto basis
EMSl = EEtml + nfEEtml;

%Energy balance by mode
EBMl = zeros(K,tstep+1);
for i=1:K
    EBMl(i,:) = sum(EMSl((i-1)*(N+1)+1:i*(N+1),:),1);
end

% dtpl = floor(tstep/40);
% XPL = XCoor;
% XPL = XPL(:,1:dtpl:tstep);
% YPL = Times';
% YPL = YPL(:,1:dtpl:tstep);
% 
% figure(5)
% 
% subplot(1,3,1)
% ZPL = EEtm(:,1:tstep);
% ZPL = ZPL(:,1:dtpl:tstep);
% surf(XPL,YPL,ZPL)
% xlabel('x')
% ylabel('t')
% title('Actual energy of the system')
% % axis([xL xR 0 FinalTime]);
% view(0,90);
% xticks(xL:(xR-xL)/K:xR)
% % shading interp
% colorbar
% 
% subplot(1,3,2)
% ZPL = nfEEtm(:,1:tstep);
% ZPL = ZPL(:,1:dtpl:tstep);
% surf(XPL,YPL,ZPL)
% xlabel('x')
% ylabel('t')
% title('Advection fluxes')
% % axis([xL xR 0 FinalTime 0 1]);
% view(0,90);
% xticks(xL:(xR-xL)/K:xR)
% % shading interp
% colorbar
% % caxis([0 5])
% 
% subplot(1,3,3)
% ZPL = EMS(:,1:tstep);
% ZPL = ZPL(:,1:dtpl:tstep);
% surf(XPL,YPL,ZPL)
% xlabel('x')
% ylabel('t')
% title('Sum of energy per Mode')
% % axis([xL xR 0 FinalTime 0 1]);
% view(0,90);
% xticks(xL:(xR-xL)/K:xR)
% % shading interp
% colorbar
% % caxis([0 5])
% 
% % sumEMS = [sum(EMS,1);sum(EMS,1)];
% % subplot(2,3,6)
% figure(6)
% XPL = kron((1:K)',ones(1,tstep));
% XPL = XPL(:,1:dtpl:tstep);
% YPL = Times(:,1:K)';
% YPL = YPL(:,1:dtpl:tstep);
% ZPL = EBM(:,1:tstep);
% ZPL = ZPL(:,1:dtpl:tstep);
% surf(XPL,YPL,ZPL)
% xlabel('x')
% ylabel('t')
% title('Sum of energy per Element')
% % axis([xL xR 0 FinalTime 0 1]);
% view(0,90);
% % shading interp
% colorbar
% % caxis([0 5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Legendre basis
figure(7)
ap =1;
for m = 1:N+1
    %Plotting the actual energy by mode for the element ap
    subplot(2,2,1)
    plot3(T,ones(length(T),1)*(m-1),EEtm(((ap-1)*(N+1))+m,:),'Color','b')
    hold on
    
    %Plotting the energy advection fluxes by mode for the element ap
    subplot(2,2,2)
    plot3(T,ones(length(T),1)*(m-1),nfEEtm(((ap-1)*(N+1))+m,:),'Color','r')
    hold on
    
    %Plotting the sum of energies by mode for the element ap
    subplot(2,2,3)
    plot3(T,ones(length(T),1)*(m-1),EMS(((ap-1)*(N+1))+m,:),'Color','g')
    hold on
    
end

%Plotting the actual energy by mode for the element ap
subplot(2,2,1)
xlabel('Time')
ylabel('Mode')
zlabel('Energy')
set(gca, 'Ydir', 'reverse')
sumle = sum(EEtm(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
plot3(T,ones(length(T),1)*(-1),sumle,'Color','k')
title(strcat('Actual Energy by mode of ',string(ap),'th domain (Black is sum)'))

%Plotting the energy advection fluxes by mode for the element ap
subplot(2,2,2)
xlabel('Time')
ylabel('Mode')
zlabel('Energy')
set(gca, 'Ydir', 'reverse')
sumlef = sum(nfEEtm(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
plot3(T,ones(length(T),1)*(-1),sumlef,'Color','k')
title(strcat('Energy Advection Flux by mode of ',string(ap),'th domain (Black is sum)'))

%Plotting the sum of energies by mode for the element ap
subplot(2,2,3)
xlabel('Time')
ylabel('Mode')
zlabel('Energy')
set(gca, 'Ydir', 'reverse')
sumlet = sum(EMS(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
plot3(T,ones(length(T),1)*(-1),sumlet,'Color','k')
title(strcat('Sum of Energies by mode of ',string(ap),'th domain (Black is sum)'))

subplot(2,2,4)
plot ( reshape(x,1,(N+1)*Elements), reshape(u,1,(N+1)*Elements), ...
          'Markersize', 8, 'LineWidth', 2 );
title('Solution')
xlabel('X')
ylabel('U')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lobatto basis
figure(8)
ap =1;
for m = 1:N+1
    %Plotting the actual energy by mode for the element ap
    subplot(2,2,1)
    plot3(T,ones(length(T),1)*(m-1),EEtml(((ap-1)*(N+1))+m,:),'Color','b')
    hold on
    
    %Plotting the energy advection fluxes by mode for the element ap
    subplot(2,2,2)
    plot3(T,ones(length(T),1)*(m-1),nfEEtml(((ap-1)*(N+1))+m,:),'Color','r')
    hold on
    
    %Plotting the sum of energies by mode for the element ap
    subplot(2,2,3)
    plot3(T,ones(length(T),1)*(m-1),EMSl(((ap-1)*(N+1))+m,:),'Color','g')
    hold on
    
end

%Plotting the actual energy by mode for the element ap
subplot(2,2,1)
xlabel('Time')
ylabel('Mode')
zlabel('Energy')
set(gca, 'Ydir', 'reverse')
sumlo = sum(EEtml(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
plot3(T,ones(length(T),1)*(-1),sumlo,'Color','k')
title(strcat('Actual Energy by lobatto mode of ',string(ap),'th domain (Black is sum)'))

%Plotting the energy advection fluxes by mode for the element ap
subplot(2,2,2)
xlabel('Time')
ylabel('Mode')
zlabel('Energy')
set(gca, 'Ydir', 'reverse')
sumlof = sum(nfEEtml(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
plot3(T,ones(length(T),1)*(-1),sumlof,'Color','k')
title(strcat('Energy Advection Flux by lobatto mode of ',string(ap),'th domain (Black is sum)'))

%Plotting the sum of energies by mode for the element ap
subplot(2,2,3)
xlabel('Time')
ylabel('Mode')
zlabel('Energy')
set(gca, 'Ydir', 'reverse')
sumlot = sum(EMSl(((ap-1)*(N+1))+1:((ap-1)*(N+1))+N+1,:),1);
plot3(T,ones(length(T),1)*(-1),sumlot,'Color','k')
title(strcat('Sum of Energies by lobatto mode of ',string(ap),'th domain (Black is sum)'))

subplot(2,2,4)
plot ( reshape(x,1,(N+1)*Elements), reshape(u,1,(N+1)*Elements), ...
          'Markersize', 8, 'LineWidth', 2 );
title('Solution')
xlabel('X')
ylabel('U')


norm(sumle,'Inf')

