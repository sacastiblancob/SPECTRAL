function [F,filterdiag] = Filter1D(N,Nc,s,V,W,alp)
% function [F] = Filter1D(N,Nc,s)
% Purpose : Initialize 1D filter matrix of size N.
% Order of exponential filter is (even) s with cutoff at Nc;
% N = Max. polynomial degree (order of interpolation)
% Nc = Number of low frequencies to remain intact
% s = Filter order
% tf = type of filter, =0 Non-Unitary (Classical), =1 Unitary (Sumedh)
% alp = alpha constant of exponential filter
    
    filterdiag = ones(N+1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Theory
    % F = VLM
    % V = vandermonde matrix
    % L = diagonal matrix with filter function
    % M = inverse of vanermonde matrix = inv(V'WV)*V'W
    % C = (V'WV)^(-1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    C = ones(N+1,1);
    C(N+1) = N/(2*N + 1);
    C = spdiags(C,0,N+1,N+1);
    
    % Initialize filter function
    
    % Exponential filter
    for i=Nc:N
        filterdiag(i+1) = exp(-alp*((i-Nc)/(N-Nc))^s);
    end
    
    F = V*spdiags(filterdiag,0,N+1,N+1)*(C*V'*W);
         
return;
end