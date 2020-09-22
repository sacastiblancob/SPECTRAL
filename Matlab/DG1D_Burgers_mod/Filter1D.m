function [F,filterdiag] = Filter1D(N,Nc,s,t,V,W,r)
% function [F] = Filter1D(N,Nc,s)
% Purpose : Initialize 1D filter matrix of size N.
% Order of exponential filter is (even) s with cutoff at Nc;
% N = Max. polynomial degree (order of interpolation)
% Nc = Number of low frequencies to remain intact
% s = Filter order
% t = type of filter, =0 Non-Unitary (Classical), =1 Unitary (Sumedh)
    
   filterdiag = ones(N+1,1);
    alp = -log(eps);
    
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
    for i=Nc:N
        filterdiag(i+1) = exp(-alp*((i-Nc)/(N-Nc))^s);
    end
    
    %t==0 means non-unitary filter
    %t==1 means Sumedh Unitary filter
    %t==2 means Lobatto Basis filter
    if t==0    
        F = V*spdiags(filterdiag,0,N+1,N+1)*(C*V'*W);
    elseif t==1
        VU = sqrt(W)*V*sqrt(C);
        F = VU*spdiags(filterdiag,0,N+1,N+1)*VU';
    elseif t==2
        VL = Basis2(r,N);
        F = VL*spdiags(filterdiag,0,N+1,N+1)*inv(VL);
    end
        
return;