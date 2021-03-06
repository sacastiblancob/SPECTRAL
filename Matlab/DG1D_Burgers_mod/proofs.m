% N = 5;
% alp = 0;
% bet = 0;
% 
% %Roots and Weights
% [x,w] = RWgll(alp,bet,N);
% W = spdiags(w,0,N+1,N+1);
% 
% %Vandermonde matrix and numerical dot product
% B = Vandermondes(x,alp,bet);
% I1 = B'*W*B;
% %for Legendre case, last term of I1 is equals to 1 + (N+1)/N
% a = I1(N+1,N+1);
% b = 1/a;
% 
% %actual inv(B'*W*B) matrix, that is not the identity but almost
% C = ones(N+1,1);
% %C(N+1) = N/(2*N + 1);      %original for Legendre case
% C(N+1) = N/(2*(N+alp)+1);
% C = spdiags(C,0,N+1,N+1);
% 
% %Unitary Sumedh version
% U = sqrt(W)*B;
% I = U*U';
% U2 = sqrt(W)*B*sqrt(C);
% I2 = U2*U2';
% % error = zeros(size((-0.4:0.05:1)'));
% % error2 = zeros(size((-0.4:0.05:1)'));
% % j=1;
% % 
% % for i=-0.4:0.05:1
% %     
% %     k=i;
% %     %Roots and Weights
% %     [x,w] = RWgll(i,k,N);
% %     W = spdiags(w,0,N+1,N+1);
% % 
% %     %Vandermonde matrix and numerical dot product
% %     B = Vandermondes(x,i,k);
% %     I1 = B'*W*B;
% %     %for Legendre case, last term of I1 is equals to 1 + (N+1)/N
% %     a = I1(N+1,N+1);
% %     b = a - 1;
% %     
% %     error(j) = a;
% %     error2(j) = 1 + (N+1+2*i)/N; %= (N+N+1+2*i)/N = (2(N+i)+1)/N
% %     j=j+1;
% % 
% % end
% % 
% % plot(-0.4:0.05:1,error)
% % hold on
% % plot(-0.4:0.05:1,error2)
% 
% %Matrix-Unitary base vs Classical Legendre Base
% 
% xp = (-1:0.01:1)';
% Bc = legendre(xp,N);
% %Bu = sqrt(w').*Bc;
% wpol = polyfit(x,sqrt(w),((N-mod(N,2))+4));
% %swint = spline(x,sqrt(w),xp);
% swint = polyval(wpol,xp);
% 
% % plot(x,sqrt(w))
% % hold on
% %plot(xp,swint)
% 
% Bu = swint.*Bc;
% 
% % % tiledlayout(2,1)
% % 
% % % nexttile
% % plot(xp,Bc)
% % title('Non-Normalized Legendre Polynomials and GLL nodes')
% % xlabel('\zeta')
% % ylabel('L(\zeta)')
% % ylim([-1.2,1.2])
% % 
% % hold on
% % scatter(x,zeros(N+1,1))
% % 
% % Bcad = legendre(xp,N+1);
% % hold on
% % plot(xp,Bcad(:,N+2),'LineStyle','--')
% % 
% % % 
% % % nexttile
% % % plot(xp,Bu)
% % % title('Matrix-Unitary Legendre')
% % % xlabel('x')
% % % ylabel('sqrt(W)*L(x)')
% % % ylim([-1.2,1.2])
% 
% %%
% %proof over some function
% 
% %the function
% %u = -sin(pi*x);
% u=ones(N+1,1);
% 
% %Non-Unitary version
% Mn = C*B'*W;
% un = Mn*u;
% unn = B*un;
% 
% %Unitary Version
% Bu = sqrt(W)*B*sqrt(C);
% Mu = Bu';
% uu = (sqrt(W)*B*sqrt(C))'*u;
% unu = (sqrt(W)*B*sqrt(C))*uu;
% 
% %Lobatto basis version
% Bl = Basis2(x,N);
% ul = inv(Bl)*u;
% unl = Bl*ul;
% 
% %Filter
% % filterdiag = ones(N+1,1);
% % alp = -log(eps);
% % s = 6;      %filter order
% % Nc = 3;     %upper cutoff
% % 
% % for i=Nc:N
% %     filterdiag(i+1) = exp(-alp*((i-Nc)/(N-Nc))^s);
% % end
% 
% %Non-Unitary Filter Matrix
% 
% %Filter with only high frequency filtered
% filterdiag = ones(N+1,1);
% filterdiag(N+1) = 0.0;
% filterdiag(N) = 0.0;
% 
% %Filter matrices
% %L = spdiags(filterdiag,0,N+1,N+1);
% L = diag(filterdiag);
% Fn = B*L*Mn;
% Fu = Bu*L*Mu;
% Fl = Bl*L*inv(Bl);
% 
% %Filtering
% Unf = Fn*u;
% Uuf = Fu*u;
% Ulf = Fl*u;
% 
% 
% % plot(x,u)
% % %legend('Original U')
% % hold on
% % plot(x,Unf)
% % %legend('U filter N-U')
% % hold on
% % plot(x,Uuf)
% % xlabel('x')
% % ylabel('u')
% % title('Filtering u=x, preserving the 2 lowest frequencies')
% % legend('U original','U filter N-U','U filter U')
% 
% % un'
% % uu'
% % ul'

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Chebyshev
% % 
% % N = 5;
% % alpc = -0.5;
% % betc = -0.5;
% % 
% % %Roots and Weights
% % [xc,wc] = RWgll(alpc,betc,N);
% % Wc = spdiags(wc,0,N+1,N+1);
% % 
% % %Chebyshev vandermonde matrix
% % xaux = (-1:0.01:1)';
% % Taux = chebyshev(xaux,N);
% % T = chebyshev(xc,N);
% % 
% % % plot(xaux,Taux)
% % % title('Non-Normalized Chebyshev Polynomials and CGL nodes')
% % % xlabel('\zeta')
% % % ylabel('T(\zeta)')
% % % ylim([-1.2,1.2])
% % % 
% % % hold on
% % % scatter(xc,zeros(N+1,1))
% % 
% % %integral proof
% % int = sum(wc.*xc.^2);  %equals to the integral between -1 and 1 of (x^2)/(sqrt(1-x^2))
% % 
% % %proofs about unitariry
% % I1 = T*T';          %Not identity but diagonal
% % T2 = sqrt(Wc)*T;    %Unitary indeed
% % I2 = T2*T2';        %Identity indeed
% % I3 = T'*Wc*T;       %Identity also!!
% % IT = T'*Wc;         %Inverse of T but without computing and inverse
% % 
% % %proofs with a function
% % y = ones(N+1,1);             %The funciton
% % umc = IT*y;          %Modal values with chevyshev
% % unc = T*umc;       %Going back to nodal, just for proof
% % 
% % umcu = T2'*y;       %Modal values for base T2
% % 
% % %modifiying Wc for force linear behaviour of W
% % Wc2 = Wc;
% % Wc2(1,1) = pi/N;
% % Wc2(N+1,N+1) = pi/N;
% % 
% % %computing new base with Wc modified
% % T3 = sqrt(Wc2)*T;   %the base
% % I4 = T3*T3';        %unitarity proof
% % umcum = inv(T3)*y;
% % 
% % 
% % %%%%%%%%%%%%%%%%%
% % %Filter
% % %%%%%%%%%%%%%%%%%
% % 
% % % filterdiag = ones(N+1,1);
% % % alp = -log(eps);
% % % s = 6;      %filter order
% % % Nc = 3;     %upper cutoff
% % % 
% % % for i=Nc:N
% % %     filterdiag(i+1) = exp(-alp*((i-Nc)/(N-Nc))^s);
% % % end
% % 
% % %Non-Unitary Filter Matrix
% % 
% % %Filter with only high frequency filtered
% % filterdiag = ones(N+1,1);
% % filterdiag(N+1) = 0.0;
% % filterdiag(N) = 0.0;
% % 
% % %Filter matrices
% % %L = spdiags(filterdiag,0,N+1,N+1);
% % L = diag(filterdiag);
% % F1 = T*L*IT;
% % F2 = T2*L*T2';
% % F3 = T3*L*inv(T3);
% % 
% % F4 = Bl*L*inv(Bl);
% % F5 = B*L*inv(B);
% % 
% % [S1,Q1,Z1] = svd(F1);
% % [S2,Q2,Z2] = svd(F2);
% % [S3,Q3,Z3] = svd(F3);
% % [S4,Q4,Z4] = svd(F4);
% % [S5,Q5,Z5] = svd(F5);


% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Lobato Basis
% 
% % Number of nodes, alpha and beta
% N = 9;
% alpha = 0.0;
% beta = 0.0;
% 
% % Computing roots and weights
% [x,w] = RWgll(alpha,beta,N);
% 
% % Computing vandermonde matrix
% [VL,M,S,K2,K4,D] = Basis2(x,N);
% W = spdiags(w,0,N+1,N+1);
% 
% % %Lobatto basis vandermonde matrix for plotting purposes
% % xaux = (-1:0.01:1)';
% % Vaux = Basis2(xaux,N);
% % 
% % plot(xaux,Vaux)
% % title('Non-Normalized Lobato Basis and GLL nodes')
% % xlabel('\zeta')
% % ylabel('T(\zeta)')
% % ylim([-1.2,1.2])
% % 
% % hold on
% % scatter(x,zeros(N+1,1))
% 
% %integral proof
% int = sum(w.*x.^2);  %equals to the integral between -1 and 1 of (x^2)/(sqrt(1-x^2))
% 
% %proofs about unitariry
% I1 = VL*VL';          %Not identity at all
% V2 = sqrt(W)*VL;     %Not unitary
% I2 = V2*V2';        %Not identity
% I3 = VL'*W*VL;        %Not identity
% IV = VL'*W;          %Not Inverse of V
% 
% %proofs with a function
% % y = ones(N+1,1);             %The funciton
% y = (sin(pi*x)+1);
% umc = inv(VL)*y;              %Modal values with Lobatto Basis
% unc = VL*umc;                 %Going back to nodal, just for proof
% 
% %Lobatto basis decomposition of y
% lode = VL.*umc'; %lode is the respective part of y stored in every lobatto basis function (by column)
% 
% % That's it: column 1 = part of y in lobatto function 1; column 2 = part of y in lobatto function 2... etc
% 
% % %Plotting the lobatto decomposition
% % plot(x,lode)
% % hold on
% % plot(x,sum(lode,2))
% 
% %Integrating Lobatto basis decomposition of y
% inte = sum(w.*(lode),1);
% 
% %%%%%%%%%%%%%%%%%
% %Filter
% %%%%%%%%%%%%%%%%%
% 
% filterdiag = ones(N+1,1);
% alp = -log(eps);
% s = 1;      %filter order
% Nc = 3;     %upper cutoff
% 
% for i=Nc:N
%     filterdiag(i+1) = exp(-alp*((i-Nc)/(N-Nc))^s);
% end
% 
% %Non-Unitary Filter Matrix
% 
% %Filter with only high frequency filtered
% % filterdiag = ones(N+1,1);
% % filterdiag(N+1) = 0.0;
% % filterdiag(N) = 0.0;
% % filterdiag(N-1) = 0.0;
% % filterdiag(N-2) = 0.0;
% 
% %Filter matrices
% %L = spdiags(filterdiag,0,N+1,N+1);
% L = diag(filterdiag);
% F1 = VL*L*inv(VL);
% 
% [S1,Q1,Z1] = svd(F1);
% 
% %Filtering
% yf = F1*y;
% 
% % With Legendre
% V = Vandermondes(x,alpha,beta);
% um = inv(V)*y;              %Modal values with Lobatto Basis
% un = V*um;
% 
% FL = V*L*inv(V);
% 
% [S,Q,Z] = svd(FL);
% 
% %Filtering
% yfl = FL*y;
% 
% % plot(x,y)
% % hold on
% % plot(x,yf)
% % hold on
% % plot(x,yfl)
% % legend('original','Lobatto Filter','Legendre Filter')
% 
% umcf = inv(VL)*yf;              %Modal values with Lobatto Basis
% uncf = VL*umcf;
% 
% umf = inv(V)*yfl;              %Modal values with Lobatto Basis
% unf = V*umf;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Strong form

N = 10;
alp = 0;
bet = 0;

%Roots and Weights
[x,w] = RWgll(alp,bet,N);
W = spdiags(w,0,N+1,N+1);
r=x;

%actual inv(B'*W*B) matrix, that is not the identity but almost
C = ones(N+1,1);
%C(N+1) = N/(2*N + 1);      %original for Legendre case
C(N+1) = N/(2*(N+alp)+1);
C = spdiags(C,0,N+1,N+1);

%Vandermonde matrix and numerical dot product
V = Vandermonde1D (0,0,N,r);
%invV = inv ( V );
invV = C*V'*W;
Dr = Dmatrix1D ( N, r, V );
Dr2 = Dr*Dr;

%integral proof
% int = sum(w.*x.^2);  %equals to the integral between -1 and 1 of (x^2)/(sqrt(1-x^2))

%the function
% u = -sin(pi*x);
u = (x-1).^2;
um = invV*u;
int = sum(w.*u);


du = Dr*u;
d2u2 = Dr2*u;

%%%
% First term (Energy of solution)
f1 = w.*(u.^2);
int1 = sum(f1);
f1m = um.^2;
int1m = sum(f1m);

%%%
% Second term (Non-linear term)
f2 = w.*(Dr*(u.^3));
int2 = sum(f2);
f2m = (invB*((Dr*(u.^3)).^0.5)).^2;
int2m = sum(f2m);
f2m2 = (invB*((3*(u.^2).*(Dr*u)).^0.5)).^2;
int2m2 = sum(f2m2);
f2m3 = (invB*((sqrt(3)*u).*(Dr*u).^0.5)).^2;
int2m3 = sum(f2m3);
f2m4 = 3*(invB*(u.*((Dr*u).^0.5))).^2;
int2m4 = sum(f2m4);

% for i=1:N+1
%     a=0;
%     if mod(i,2)==0
%         a = a + um(i);
%     else
%         a = a - um(i);
%     end
%     a
% end




