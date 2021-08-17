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

% i=1;
for i=0:(ne-1)

    NM = reshape(RES(ind+i,1),[nk nn]);
    KM = reshape(RES(ind+i,2),[nk nn]);
    
    %MER with sign
    EM = reshape(RES(ind+i,7),[nk nn]);
    EMsign = EM;
    EMsign(abs(EM)<=dt) = 0;
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
    xlabel('P. Order')
    ylabel('Elements')
    zlabel('MER')
    view(0,90)
    axis([9 48 1 16 -0.5 0.5])
%     c=hot(10);
%     colormap(c);
%     colormap parula
%     colormap winter
    mycolors = [51/255 255/255 51/255;0/255 255/255 255/255; 0/255 0/255 255/255];
    colormap(mycolors);
    shading interp
    caxis([-0.03 0.03])
    title(strcat('\nu = ',string(epsilons(i+1))))

end











