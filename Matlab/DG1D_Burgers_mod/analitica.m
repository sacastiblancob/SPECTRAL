%Solver of the analitic solution of the Burgers Ecuation, using
%Gauss-Legendre numerical integration for 11 nodes.
%Notation:          analitica(x,t,nu)
%Where: x=vector with x values ; t=scalar with time step;
%       nu=kinematic viscosity
%

function I=analitica(x,t,nu)

[eta,w] = JacobiGQ(0,0,100);

a = -1;
b = -a;
%I=zeros(1,max(size(x)));
s = size(x);
I = zeros(s(1),s(2));
for i=1:s(1)*s(2)    %i=1:max(size(x))
    up=((b-a)/2)*sum(sin(pi.*(x(i)-(((b-a)*eta+(b+a))/2))).*exp(-cos(pi*(x(i)-(((b-a)*eta+(b+a))/2)))/(2*pi*nu)).*exp(-((((b-a)*eta+(b+a))/2).^2)/(4*nu*t)).*w);
    down=((b-a)/2)*sum(exp(-cos(pi*(x(i)-(((b-a)*eta+(b+a))/2)))/(2*pi*nu)).*exp(-((((b-a)*eta+(b+a))/2).^2)/(4*nu*t)).*w);
    I(i)=-(up/down);
end

end





