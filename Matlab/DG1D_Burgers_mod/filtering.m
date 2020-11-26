function uf = filtering(u,F)
%Function for apply a filter over SEM non-linear solution

% uvf = uv;
% for k = 1:K
%     uvf((Np-1)*(k-1)+1:(Np-1)*(k)+1) = F*uv((Np-1)*(k-1)+1:(Np-1)*(k)+1);
% end
%u = vec2matrix(uv,K);
siz = size(u);
K = siz(2);
Np = siz(1);

%Applying filtering
uf = F*u;

% %Ensure continuity at element interfaces
% for k = 2:K
%     uf(Np,k-1) = (uf(1,k)/2) + (uf(Np,k-1)/2);
%     uf(1,k) = uf(Np,k-1);
%     %uf(Np,k-1) = u(Np,k-1);        %uncomment for not filter the boundaries of the elements
%     %uf(1,k) = u(1,k);              %uncomment for not filter the boundaries of the elements
% end

end