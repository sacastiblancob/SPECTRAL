function v = mat2vec(M)
%This function takes matrix (as DG-SEM matrices) and return the vectors for
%work with continuous SEM methods

siz = size(M);
K = siz(2);
Np = siz(1);

Ns=((K-1)*(Np-1))+Np;
v = zeros(Ns,1);
v(1) = M(1,1);
v(Ns) = M(Np,K);
for k = 1:K
    v((Np-1)*(k-1)+2:(Np-1)*(k)+1) = M(2:Np,k);
end

end