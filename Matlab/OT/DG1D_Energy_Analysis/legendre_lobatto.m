function L = legendre_lobatto(x,K)
%% Return the value of the Legendre polynomials of order <=K at the points x
%Nidia's work
%Modification for could evaluate polynomials over any vector x
%computes the orthogonal (not orthonormal) version of the Legendre
%polynomials

Nx=length(x);
%Nx = K+1;
L=zeros(Nx,K+1);

L(:,1)=ones(Nx,1);
L(:,2)=x;

for ii=1:K-1
    L(:,ii+2) = ((2*ii+1)/(ii+1)).*x.*L(:,ii+1) - ((ii)/(ii+1)).*L(:,ii);
end

% for j=1:K+1
%     L(:,j)=L(:,j)/sqrt(2.0/(2.0*(j-1)+1));
% end

end
