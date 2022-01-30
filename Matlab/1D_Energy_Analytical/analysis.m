%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading results
% load('RES_ana.mat')
load('RES_ana_full.mat')

%variables
% Nes = 9:3:50;
Nes = 10:50;
% Kes = [1 2 3 4 7 10 13 16];
Kes = 1:16;
epsilons = [0.05 0.01 0.0050 0.001 0.0005 0.00025];

dt = 0.005;

%number of variables
nn = length(Nes);
nk = length(Kes);
ne = length(epsilons);

% N vs Elements vs MAE
ind = 1:ne:nn*nk*ne;

dts = zeros(size(RES(:,1)));
for i=1:length(dts)
    k = RES(i,1);
    N = RES(i,2);
    epsilon = RES(i,3);
    
    % Generate simple mesh
    xL = -1;
    xR = 1;
    [Nv, VX, K, EToV] = MeshGen1D(xL,xR,k);
    
    %
    %  Build coordinates of all the nodes
    %
    r = JacobiGL ( 0, 0, N );
    va = EToV(:,1)'; 
    vb = EToV(:,2)';
    %normal version
    x = ones(N+1,1) * VX(va) + 0.5 * (r+1) * (VX(vb)-VX(va));
    
    xmin = min(abs(x(1,:)-x(2,:)));
    CFL = 1.0;
    umax = 1.0;
    dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));
    
    dts(i) = dt;
    
end

% i=1;
for i=0:(ne-1)

    NM = reshape(RES(ind+i,1),[nk nn]);
    KM = reshape(RES(ind+i,2),[nk nn]);
    
    %MER with sign
    EM = reshape(RES(ind+i,7),[nk nn]);
    DTM = reshape(dts(ind+i),[nk nn]);
    EMsign = sign(EM);
    EMsign(abs(EM)<= 0.05) = 0;
%     EMsign = 100*EMsign;
    
    %NORM inf
    ENIM = reshape(RES(ind+i,6),[nk nn]);
    ENIMlog = log(abs(ENIM));
    % EMlog = s*EMlog;

%     figure(i+1)
% 
%     subplot(1,2,1)
%     surf(NM,KM,EM)
%     xlabel('Nodes')
%     ylabel('Elements')
%     zlabel('MER')
%     axis([9 48 1 16 -0.5 0.1])
%     shading interp
% 
%     subplot(1,2,2)
%     surf(NM,KM,ENIMlog)
%     xlabel('Nodes')
%     ylabel('Elements')
%     zlabel('log - Error NormInf')
%     view(0,90)
%     axis([9 48 1 16 -8 1])
%     colorbar
%     caxis([-8 1])
% 
%     sgtitle(strcat('\nu = ',string(epsilons(i+1))))
    
    figure(10)
    
    subplot(2,3,i+1)
    surf(NM,KM,EM.*abs(EMsign))
%     surf(NM,KM,EM)
    view(0,90)
    axis([10 48 1 16 -0.5 0.5])
    xlabel('P. Order')
    ylabel('Elements')
    zlabel('MER')
%     c=hot(10);
%     colormap(c);
%     colormap bone
%     colormap winter
%     mycolors = [51/255 255/255 51/255;0/255 255/255 255/255; 0/255 0/255 255/255];
    mycolors = [150/255 150/255 150/255; 220/255 220/255 220/255; 50/255 50/255 50/255];
    colormap(mycolors);
    shading interp
    caxis([-0.01 0.01])
    ttl = title(strcat('R_e = ',string(2/epsilons(i+1))));
    
%     figure(100)
%     
%     subplot(2,3,i+1)
%     surf(NM,KM,DTM)
%     view(0,90)
% %     axis([9 48 1 16 -0.5 0.5])
%     xlabel('P. Order')
%     ylabel('Elements')
%     zlabel('MER')
% %     c=hot(10);
% %     colormap(c);
% %     colormap parula
% %     colormap winter
%     colormap(mycolors);
%     shading interp
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NUMERICAL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % loading results
% % load('RES_ana.mat')
% load('RES_num_save46_2.mat')
% 
% % %variables
% % % Nes = 9:3:50;
% % Nes = 10:50;
% % % Kes = [1 2 3 4 7 10 13 16];
% % Kes = 1:16;
% % epsilons = [0.05 0.01 0.0050 0.001 0.0005 0.00025];
% 
% dt = 0.005;
% 
% % %number of variables
% % nn = length(Nes);
% % nk = length(Kes);
% % ne = length(epsilons);
% 
% % N vs Elements vs MAE
% ind = 1:ne:nn*nk*ne;
% 
% % i=1;
% for i=0:(ne-1)
% 
%     NM = reshape(RES(ind+i,1),[nk nn]);
%     KM = reshape(RES(ind+i,2),[nk nn]);
%     
%     %MER with sign
%     EM = reshape(RES(ind+i,7),[nk nn]);
%     EMsign = sign(EM);
%     EMsign(abs(EM)<=0.05) = 0;
% %     EMsign = 100*EMsign;
%     
%     %NORM inf
%     ENIM = reshape(RES(ind+i,6),[nk nn]);
%     ENIMlog = log(abs(ENIM));
%     % EMlog = s*EMlog;
% 
% %     figure(i+1)
% % 
% %     subplot(1,2,1)
% %     surf(NM,KM,EM)
% %     xlabel('Nodes')
% %     ylabel('Elements')
% %     zlabel('MER')
% %     axis([9 48 1 16 -0.5 0.1])
% %     shading interp
% % 
% %     subplot(1,2,2)
% %     surf(NM,KM,ENIMlog)
% %     xlabel('Nodes')
% %     ylabel('Elements')
% %     zlabel('log - Error NormInf')
% %     view(0,90)
% %     axis([9 48 1 16 -8 1])
% %     colorbar
% %     caxis([-8 1])
% % 
% %     sgtitle(strcat('\nu = ',string(epsilons(i+1))))
%     
%     figure(11)
%     
%     subplot(2,3,i+1)
% %     EM(isnan(EM))=99;
%     EMsign(isnan(EMsign))=1;
%     surf(NM,KM,EM.*abs(EMsign))
% %     EM(isnan(EM))=99;
% %     contour(NM,KM,EM.*abs(EMsign),100)
%     view(0,90)
%     xlh = xlabel('P. Order');
%     ylh = ylabel('Elements');
%     zlh = zlabel('MER');
%     ttl = title(strcat('\nu = ',string(epsilons(i+1))));
% %     xlh.Position = [28.5 -0.7986 -0.5000];
% %     ylh.Position = [2.6279 7.5 -0.5];
% %     zlh.Position = [9 1 1];
% %     ttl.Position = [28.5 16.2435 0];
%     xlh.Units =  'normal';
%     zlh.Units =  'normal';
%     ylh.Units =  'normal';
%     ttl.Units =  'normal';
%     xl = 9;
%     xr = 48;
%     yl = 1;
%     yr = 16;
%     zl = -100;
%     zr = 100;
%     axis([xl xr yl yr zl zr])    
% %     c=hot(10);
% %     colormap(c);
% %     colormap parula
% %     colormap winter
%     mycolors = [51/255 255/255 51/255;0/255 255/255 255/255; 0/255 0/255 255/255];
%     colormap(mycolors);
%     shading interp
%     caxis([-0.03 0.03])
%     
% 
% end











