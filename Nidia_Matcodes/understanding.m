clear all
clc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROOFS ABOUT MODAL BASIS FUNCTIONS
% *Modal boundary adapted basis
% Sergio Castiblanco
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%entries
n = 25;                          %number of nodes
p = n-1;                        %max polynomial order
alp = 0;                        %jacobi alpha, 0.0 for Legendre
bet = 0;                        %jacobi beta, 0.0 for Legendre

%Legendre
%x = JacobiGL(alp,bet,p);    %GLL nodes and weights
[x,w] = gll(p);
B = Vandermondes(x,alp,bet);    %Legendre Vandermonde matrix
W = diag(w);
RW = diag(sqrt(w));

%modal basis
[V,M,S,K2,K4,D]=Basis2(x,p);                %modal base
[VN,MN,SN,K2N,K4N,DN] = Basis2norm(x,p);    %normalized modal base

%filter matrix
eps = -log(1.0e-16);            %Machine precision according to Diamesis 2008
pf = 30;                        %filter order
xauxf = 1:n;
xauxf = ((xauxf - 1.0) / (n - 1.0));    %auxiliar x for compute filter
lf = exp(-eps * xauxf.^pf);             %vector with filter values
LF = diag(lf);                          %diag. matrix for filter

% for i=1:n
%     plot(x,V(:,i))
%     drawnow
%     hold on
% end

%filtering matrices
F1 = B*LF*inv(B);
F2 = V*LF*inv(V);

%unitarity
U1 = B'*W*B;
U2 = V'*W*B;
U3 = B'*W*V;
U4 = V'*W*V;

%norms = [2/3 2/3 2/5 2/(5+16) 2/(5+16+24) 2/(5+16+24+32) 2/(5+16+24+32+40) 2/(5+16+24+32+40+48) 2/(5+16+24+32+40+48+56)];
%[0,2,5,9,14]
%eigen-analisys
%[V,D] = eig(U4);

%SVD
[U1,S1,Z1] = svd(F1);
[U2,S2,Z2] = svd(F2);


kas = 15;
times = 20;
norm1 = zeros(kas,times+1);
norm1(:,1) = 1;
norm2 = zeros(kas,times+1);
norm2(:,1) = 1;

for k=1:kas
    
    yi = sin(2*k*x*pi);
    yo1 = yi;
    yo2 = yi;
       
    for j=1:times
        yf1 = F1*yo1;
        yf2 = F2*yo2;
        yn1 = norm(yf1)/norm(yi);
        yn2 = norm(yf2)/norm(yi);
        yo1 = yf1;
        yo2 = yf2;
        norm1(k,j+1) = yn1;
        norm2(k,j+1) = yn2;
    end
    
    subplot(1,2,1)
    plot(0:times,norm1(k,:))
    title('Legendre Basis ; u0=sin(2*k*pi*x)')
    xlabel('Times of filtering')
    ylabel('Relative Norm')
    drawnow
    hold on
    
    subplot(1,2,2)
    title('Modal Basis ; u0=sin(2*k*pi*x)')
    plot(0:times,norm2(k,:))
    xlabel('Times of filtering')
    ylabel('Relative Norm')
    drawnow
    hold on
       
end





