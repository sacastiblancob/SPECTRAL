open('normuv.mat')
normuv = ans.normuv;

figure(1)
semilogy(normuv(:,1),normuv(:,3),'k-','LineWidth',1)
hold on
semilogy(normuv(:,10),normuv(:,12),'k--','LineWidth',1)
hold on
semilogy(normuv(:,4),normuv(:,6),'k:','LineWidth',1.5)
hold on
semilogy(normuv(:,13),normuv(:,15),'k-.','LineWidth',1)
hold on
legend('Classical Approach \nu=0.0005','Classical Approach \nu=0.00025','Energy-Based Filter\nu=0.0005','Energy-Based Filter\nu=0.00025','Location','Southeast')
grid on
xlabel('Time (s)')
ylabel('Absolute Inf Error')

open('solut.mat')
solut = ans.solut;

figure(2)

% plot(solut{1,1},solut{1,2},'k-.')
% hold on
plot(solut{1,3},solut{1,4},'k:')
hold on
plot(solut{1,5},solut{1,6},'k--')
hold on
plot(solut{1,7},solut{1,8},'k-')
grid on







