function uvf = filtering(uv,F,K,Np)
%Function for apply a filter over SEM non-linear solution

% uvf = uv;
% for k = 1:K
%     uvf((Np-1)*(k-1)+1:(Np-1)*(k)+1) = F*uv((Np-1)*(k-1)+1:(Np-1)*(k)+1);
% end
u = vec2matrix(uv,K);
uf = F*u;
for k = 2:K
    uf(Np,k-1) = (uf(1,k)/2) + (uf(Np,k-1)/2);
    %uf(Np,k-1) = u(Np,k-1);        %uncomment for not filter the boundaries of the elements
    %uf(1,k) = u(1,k);              %uncomment for not filter the boundaries of the elements
end
uvf = mat2vec(uf);

% Ns = length(uv);
% ap = 0;
% UN = zeros(Ns,1);
% UT = zeros(Np,1);
% for i=1:K
% 	for j=1:Np
%         UT(j)=uv(Np*(i-1)+j-ap);
%     end
% 	UF=F*UT;
%     if (i ~= K)
%         UF(Np)=UF(Np)/2.0;
%     end
% 	if (i ~= 1)
% 	    UF(1)=UF(1)/2.0;
%     end
% 	for j=1:Np
%         UN(Np*(i-1)+j-ap)=UN(Np*(i-1)+j-ap)+UF(j);
%     end
% 	ap=ap+1;
% end
% uvf2=UN;

end