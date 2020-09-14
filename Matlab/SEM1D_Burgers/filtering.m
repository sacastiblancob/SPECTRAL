function uvf = filtering(uv,F,K,Np)
%Function for apply a filter over SEM non-linear solution

uvf = uv;
for k = 1:K
    uvf((Np-1)*(k-1)+1:(Np-1)*(k)+1) = F*uv((Np-1)*(k-1)+1:(Np-1)*(k)+1);
end

end