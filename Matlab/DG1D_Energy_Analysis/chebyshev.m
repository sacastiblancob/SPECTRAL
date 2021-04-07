function T = chebyshev(x,N)
%Returns de Vandermonde matrix T for the chebyshev polynomials at points x
%with a basis of degree up to N

Nx=length(x);
T=zeros(Nx,N+1);

T(:,1)=ones(Nx,1);
T(:,2)=x;

for ii=1:N-1
    T(:,ii+2) = 2.*x.*T(:,ii+1) - T(:,ii);
end

% for j=1:N+1
%     if (j==1 || j==N+1)
%         T(:,j) = T(:,j)/sqrt(pi);
%     else
%         T(:,j) = T(:,j)/sqrt(pi/2);
%     end
% end

end