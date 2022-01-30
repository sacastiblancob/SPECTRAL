%Solver of the analitic solution of the Burgers Ecuation, using
%Gauss-Legendre numerical integration for 11 nodes.
%Notation:          analitica(x,t,nu,eta,w,a)
%Where: x=vector with x values ; t=scalar with time step;
%       nu=kinematic viscosity
% for t>0

function I=analitica(x,t,nu,etai,wi,ai)

if nu~= 0
    bi = -ai;
    %I=zeros(1,max(size(x)));
    s = size(x);
    I = zeros(s(1),s(2));

    for i=1:s(1)*s(2)    %i=1:max(size(x))
        h = (((bi-ai)*etai+(bi+ai))/2);
        xmh = (x(i)-h);
    %     f1 = exp(-1/(2*pi*nu) + Initial_condition_int(xmh,nu)/(-2*pi*nu));
        f = exp(Initial_condition_int(xmh,nu)/(-2*pi*nu));
        e = exp(-(h.^2)/(4*nu*t));
        ef = e.*f;
    %     ef1 = e.*f1;
        up = ((bi-ai)/2)*sum(ef.*Initial_condition(xmh,nu).*wi);
        down = ((bi-ai)/2)*sum(ef.*wi);
        I(i) = (up/down);

    %      up=((b-a)/2)*sum(sin(pi.*(x(i)-(((b-a)*eta+(b+a))/2))).*exp(-cos(pi*(x(i)-(((b-a)*eta+(b+a))/2)))/(2*pi*nu)).*exp(-((((b-a)*eta+(b+a))/2).^2)/(4*nu*t)).*w);
    %     down=((b-a)/2)*sum(exp(-cos(pi*(x(i)-(((b-a)*eta+(b+a))/2)))/(2*pi*nu)).*exp(-((((b-a)*eta+(b+a))/2).^2)/(4*nu*t)).*w);
    %     I(i)=-(up/down);
    end
else
    s = size(x);
    I = zeros(s(1),s(2));
    for i=1:s(1)*s(2)    %i=1:max(size(x))
        xaux = x(i);
        minres = 1E-14;
%         x0 = xaux;
        x0 = sign(xaux);
        for j=1:100
          f0 = x0 + Initial_condition(x0)*t - xaux; %Calculating the value of function at x0
          f0_der = 1 + Initial_condition_prima(x0)*t; %Calculating the value of function derivative at x0
          y = x0 - f0/f0_der; % The Formula
          err=abs(y-x0);
        if err<minres %checking the amount of error at each iteration
            break
        end
          x0=y;
        end
        I(i) = Initial_condition(y);
    end
    
end
end





