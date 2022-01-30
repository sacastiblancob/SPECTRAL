%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sergio Castiblanco
% Pontificia Unviersidad Javeriana
% 2021
%
% Solver of 1D-Advection-Diffusion-Reaction Equation with
% Fourier Galerking Methodology
%    du/dt + a du/dx = nu d2u/dx2 + re u
%

%%%%%%%%
%% Initial Configuration
% Number of grid points
N = 101;

% Boundaries
xL = -pi;
xR = pi;

% Constants
a = 0.5;
nu = 0.0;
re = 0.0;

%% Grid Construction
% Space
dx = (xR-xL)/N;
x = -(xR-xL)/2:dx:(xR-xL)/2-dx;
x = x';
if mod(N,2) == 0
    K = (2*pi/(xR-xL))*(-N/2:N/2-1);
    K = K';
    K = fftshift(K);
else
    K = (2*pi/(xR-xL))*floor(-N/2+1:N/2);
    K = K';
    K = ifftshift(K);
end


% Time grid
to = 0.0;
tf = pi;
dt = 0.01;
t = to:dt:tf;
Ts = length(t);

%% Initial Condition
uo = sin(2*(2*pi/(xR-xL))*x);
% uo = (1/sqrt(4*pi*nu*to))*exp(-((x+a).^2)/(4*nu*to));

umo = fft(uo);

%% SOLVING
% Integrating factor
gamma = -i*K*a - nu*K.^2 - re;
I = exp(gamma*dt);

for time = 1:Ts
    um = I.*umo;
        
    u = ifft(um);
    ua = sin(2*(2*pi/(xR-xL))*(x - a*time*dt));
    
    umo = um;
    
%     if a>=0
%         u(1) = u(end);
%     else
%         u(end) = u(1);
%     end
    
    plot(x,u)
    hold all
    plot(x,ua)
    xlim([xL xR])
    ylim([-1.1 1.1])
    hold off
    drawnow
    
%      pause(10*dt)
    pause(0.01)
    
    
end






