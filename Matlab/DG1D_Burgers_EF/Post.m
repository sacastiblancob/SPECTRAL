% %
% % Post-Processing
% %

% figure(1)
%         
% subplot(1,2,1)
% plot(xv,uv_ebf);
% hold all
% plot(xv,uv_svv);
% plot(xv,uv_cf);
% plot(xv,uva)
% xlim([xL xR]);
% ylim([ulbo-0.1 ulup+0.1]);
% title('Solution u, ua')
% xlabel('X')
% ylabel('U')
% legend('EBF','SVV','CF','Analytic')
% 
% subplot(1,2,2)
% semilogy(xv,Error_ebf);
% hold all
% semilogy(xv,Error_svv);
% semilogy(xv,Error_cf);
% xlim([xL xR]);
% ylim([1E-16 1]);
% title('Error Structure |ua - u|')
% xlabel('X')
% ylabel('log(Error)')
% legend('EBF','SVV','CF')
% 
% Errorv_exp = reshape(Error_exp,1,[]);
% % Errorv_van = reshape(Error_van,1,[]);
% Errorv_cf = reshape(Error_cf,1,[]);
% Errorv_svv = reshape(Error_svv,1,[]);
% 
% figure(8)
% semilogy(xv,Errorv_exp)
% hold all
% % semilogy(xv,Errorv_van)
% semilogy(xv,Errorv_cf)
% semilogy(xv,Errorv_svv)
% legend('EBF EXP','CF','SVV')
% title('Error')
% xlabel('X')
% ylabel('|ua-u|')
% grid on
% 
% figure(9)
% plot(T,Err2_t_exp)
% hold all
% % plot(T,Err2_t_van)
% plot(T,Err2_t_cf)
% plot(T,Err2_t_svv)
% legend('EBF EXP','CF','SVV')
% 
% figure(10)
% plot(T,Errinf_t_exp)
% hold all
% % plot(T,Errinf_t_van)
% plot(T,Errinf_t_cf)
% plot(T,Errinf_t_svv)
% legend('EBF EXP','CF','SVV')
% title('Inf Error')
% xlabel('X')
% ylabel('Norm-Inf|ua-u|')
% grid on
% 
% figure(101)
% semilogy(T,Errinf_t_exp)
% hold all
% % plot(T,Errinf_t_van)
% semilogy(T,Errinf_t_cf)
% semilogy(T,Errinf_t_svv)
% legend('EBF EXP','CF','SVV')
% title('Inf Error')
% xlabel('X')
% ylabel('Norm-Inf|ua-u|')
% grid on
% 
% figure(11)
% plot(T,Errord_t_exp)
% hold all
% % plot(T,Errinf_t_van)
% plot(T,Errord_t_cf)
% plot(T,Errord_t_svv)
% legend('EBF EXP','CF','SVV')
% title('Mean order of the solution')
% xlabel('Time')
% ylabel('Order')
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Error_Res.mat

Errorv_exp = reshape(Error_exp,1,[]);
Errorv_svv = reshape(Error_svv,1,[]);
Errorv_cf_4 = reshape(Error_cf_4,1,[]);
Errorv_cf_8 = reshape(Error_cf_8,1,[]);
Errorv_cf_12 = reshape(Error_cf_12,1,[]);
Errorv_cf_16 = reshape(Error_cf_16,1,[]);
% Errorv_cf_20 = reshape(Error_cf,1,[]);

% Inf Error plot
figure(111)

subplot(1,3,1)
semilogy(T,Errinf_t_exp,'-','Color','black')
hold all
semilogy(T,Errinf_t_svv,'-.','Color','black')
semilogy(T,Errinf_t_cf_4,'--','Color','black')
semilogy(T,Errinf_t_cf_8,'--','Color',[50/255 50/255 50/255])
semilogy(T,Errinf_t_cf_12,'--','Color',[100/255 100/255 100/255])
semilogy(T,Errinf_t_cf_16,'--','Color',[150/255 150/255 150/255])
% semilogy(T,Errinf_t_cf)
legend('EBF','SVV','CF s=4','CF s=8','CF s=12','CF s=16','Location','Southeast')
title('a) Log-Inf-Error, nu=2.5E-4')
xlabel('Time')
ylabel('||u_a - u||_\infty')
grid on

subplot(1,3,2)
plot(T,Errinf_t_exp,'-','Color','black')
hold all
plot(T,Errinf_t_svv,'-.','Color','black')
plot(T,Errinf_t_cf_4,'--','Color','black')
plot(T,Errinf_t_cf_8,'--','Color',[50/255 50/255 50/255])
plot(T,Errinf_t_cf_12,'--','Color',[100/255 100/255 100/255])
plot(T,Errinf_t_cf_16,'--','Color',[150/255 150/255 150/255])
% semilogy(T,Errinf_t_cf)
% legend('EBF','SVV','CF s=4','CF s=8','CF s=12','CF s=16','Location','Northwest')
title('b) Inf-Error, nu=2.5E-4')
xlabel('Time')
ylabel('||u_a - u||_\infty')
grid on

% Order of mean error
% figure(112)
subplot(1,3,3)
plot(T,Errord_t_exp,'-','Color','black')
hold all
plot(T,Errord_t_svv,'-.','Color','black')
plot(T,Errord_t_cf_4,'--','Color','black')
plot(T,Errord_t_cf_8,'--','Color',[50/255 50/255 50/255])
plot(T,Errord_t_cf_12,'--','Color',[100/255 100/255 100/255])
plot(T,Errord_t_cf_16,'--','Color',[150/255 150/255 150/255])
% plot(T,Errord_t_cf)
% legend('EBF','SVV','CF s=4','CF s=8','CF s=12','CF s=16','Location','Southeast')
title('c) Mean order of error, nu=2.5E-4')
xlabel('Time')
ylabel('O(u_a - u)')
grid on

% Error Structure
figure(113)
semilogy(xv,Errorv_exp,'-','Color','black')
hold all
semilogy(xv,Errorv_svv,'-.','Color','black')
% semilogy(xv,Errorv_cf_4,'-','Color','black')
% semilogy(xv,Errorv_cf_8,'-','Color',[50/255 50/255 50/255])
% semilogy(xv,Errorv_cf_12,'-','Color',[100/255 100/255 100/255])
% semilogy(xv,Errorv_cf_16,'-','Color',[150/255 150/255 150/255])
semilogy(xv,Errorv_cf_16,'--','Color',[100/255 100/255 100/255])
% semilogy(xv,Errorv_cf_20)
% legend('EBF','SVV','CF s=4','CF s=8','CF s=12','CF s=16','Location','Northwest')
legend('EBF','SVV','CF s=16','Location','South','FontSize',8)
title('Error Structure, nu=2.5E-4')
% title('Error')
xlabel('x')
ylabel('|u_a - u|')
grid on



