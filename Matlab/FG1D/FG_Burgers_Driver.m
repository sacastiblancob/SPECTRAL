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
N = 51;

% Boundaries
xL = -1;
xR = 1;

% Constants
nu = 0.0;

% Filter?
cf = 0;

% Hyperviscosity?
hv = 1;

%% Grid Construction
% Space
dx = (xR-xL)/N;
x = -(xR-xL)/2:dx:(xR-xL)/2-dx;
x = x';
if mod(N,2) == 0
    K = (2*pi/(xR-xL))*(-N/2:N/2-1);
else
    K = (2*pi/(xR-xL))*(ceil(-N/2):floor(N/2));
end
K = K';
K = fftshift(K);

%% Initial Condition
u = -sin(pi*x);

%% Time
to = 0.0;
tf = 3/pi;
xmin = (xR-xL)/N;
CFL = 0.01;
umax = max(max(abs(u)));
dt = CFL* min(xmin/umax,xmin^2/sqrt(nu));
t = to:dt:tf;
Ts = length(t);

%% Filter

if cf == 1
    alp = -log(eps());
    p = 10;             %Filter order
    Kmax = max(abs(K));
    Kf = ones(N,1);
    for j=K
        Kf = exp(-alp*(K/Kmax).^(2*p));
    end
end

%% Hyperviscous operator

if hv == 1
    p = 2;
    alp = -log(eps());
end

%% SOLVING
%
%  Runge-Kutta coefficients and residual storage  
%
resu = zeros(N,1);
RK4_coef;

for time = 1:Ts
    
    for INTRK = 1:5
        um = fft(u);
        dum = i*K.*um;
        ddum = -(K.^2).*um;
        du = ifft(dum);
        ddu = ifft(ddum);
        
        if hv~=1
            rhsu = -u.*du + nu*ddu;
        else
            d2pum = (i^(2*p))*(K.^(2*p)).*um;
            rhsu = -u.*du + nu*ddu -(alp*(-1^p)/(dt*N^(2*p))).*d2pum;
        end
               
        %Solving U with RK4
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
        
        if cf==1
            um = fft(u);
            um = um.*Kf;
            u = ifft(um);
        end
       
    end
    
    plot(x,real(u))
    hold all
%     plot(x,ua)
    xlim([xL xR])
    ylim([-1.1 1.1])
    hold off
    drawnow
    
%      pause(10*dt)
%     pause(0.01)
    pause
    
    
end






