function M = vec2matrix(v,K)
%This function takes a vector (as the ones used in continuous SEM) and
%return the matrix aranged by spectral elements (as the matrices used in
%DG-SEM). K is the number of elements (number of columns of M)

Ns = length(v);
Np = floor(Ns/K)+1;
M = zeros(Np,K);
M(1) = v(1);

for k = 1:K
    M(2:Np,k) = v((Np-1)*(k-1)+2:(Np-1)*(k)+1);
end
for k = 2:K
    M(1,k) = M(Np,k-1);
end

end