function [P] = legendreM(x,N)

    %computing the size of x
    M = length(x);
    P = zeros(N,M);
    P(1,:) = 1;
    P(2,:) = x;
    for k = 1:N-2
        P(k+2,:) = ((2*k+1)/(k+1)).*x.*P(k+1,:)' - (k/(k+1)).*P(k,:)';
    end
    P = P';

return