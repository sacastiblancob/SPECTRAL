function [x,w,D,ii] = gll(Nin)
% Calculate the N GLL points and quadrature weights
% Also return the differentiation matrix
maxIter=150;
Np1=Nin+1;

x=-1.0*cos(pi*(0:Np1-1)/(Np1-1))';
%x=linspace(-1,1,Np1)'; x=sign(x).*sqrt(abs(x));

for ii=1:maxIter
   xold=x; 
   
   Lk=ones(Np1,1);
   Lkp1=xold;
   
   for kk=2:Np1-1
      Lkm1=Lk;
      Lk=Lkp1;
      Lkp1=((2*kk-1)*xold.*Lk - (kk-1)*Lkm1)/kk;
   end
   
   x = xold - (xold.*Lkp1-Lk)./(Np1*Lkp1);
   
   if(max(abs(x-xold)) < 1.1*eps) 
       break;
   end
end

if(ii == maxIter)
   error('MAXIMUM ITERATION REACHED IN GLL!!'); 
   
end

w= 2.0./((Np1-1)*Np1*Lkp1.^2)';

D=zeros(Np1,Np1);

for j=1:Np1

    xj=x(j);
    for k=[1:j-1,j+1:Np1]
       D(k,j)=(Lkp1(k)/Lkp1(j))/(x(k)-xj); 
    end
    
end

D(1,1)=-(Np1*(Np1-1))/4.0;
D(Np1,Np1)=(Np1*(Np1-1))/4.0;

end

