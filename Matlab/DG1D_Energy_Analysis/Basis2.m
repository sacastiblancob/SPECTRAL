%---------------------------------------------------------------
% Calculate the basis functions, mass, and first/second/fourth
% order stiffness matrices using numerical quadrature for the
% Legendre basis with double roots at x=-1,1.
%
% \author   Kris Rowe
% \date     Oct 2018
% \version  1.0
%
%---------------------------------------------------------------
function [V,M,S,K2,K4,D]=Basis2(x,N)

Nq=1024;     %The number of quadrature points ued to compute matrix entries
tol=1e-14;  %The threshold for deciding a matrix entry is zero
[xi,w,D]=gll(Nq);

%- Basis Functions at Quadrature Points Used to Calculate Matrix Entries
L=legendre_lobatto(xi,N+4);
V=zeros(Nq+1,N+1);
V(:,1)=(1-xi)/2;
V(:,2)=(1+xi)/2;

for ii=2:N
   V(:,ii+1)=                       L(:,ii-1) ...
                                   -L(:,ii+1) ;
   V(:,ii+1)=V(:,ii+1)/(sqrt(4*ii-2));
end

%V0=V;
%V=[V0(:,1) V0(:,3:end) V0(:,2)];

dV=D*V;     %Collocation Derivative
W=diag(w);  %Quadrature Weights

M  =V'*(W*V);      M(abs(M)<tol)  =0.0;     %M=sparse(M);
S  =V'*(W*dV);     S(abs(S)<tol)  =0.0;     %S=sparse(S);
K2  =dV'*(W*dV);   K2(abs(K2)<tol) =0.0;    %K2=sparse(K2);
K4  =speye(N+1);

%- Basis Functions at the Grid Passed Into the Function
L=legendre_lobatto(x,N+4);
V=zeros(length(x),N+1);
%[xi,w,D]=gll(N);
V(:,1)=(1-x)/2;
V(:,2)=(1+x)/2;

for ii=2:N
   V(:,ii+1)=                       L(:,ii-1) ...
                                   -L(:,ii+1) ;
   V(:,ii+1)=V(:,ii+1)/(sqrt(4*ii-2));
end

%V0=V;
%V=[V0(:,1) V0(:,3:end) V0(:,2)];

end
