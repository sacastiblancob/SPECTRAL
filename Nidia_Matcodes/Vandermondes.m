function [PL] = Vandermondes(x,alpha,beta)
    % function [P] = JacobiP(x,alpha,beta,N)
    % Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
    % (alpha+beta <> -1) at points x for order N and lower and returns
    % P[length(xp):length(xp))], vandermonde matrix
    % Note : They are normalized to be orthonormal.
    % Turn points into row if needed.
    
    N = length(x)-1;
    xp = x;
    dims = size(xp);
    if(dims(2)==1)
        xp = xp'; 
    end
    PL = zeros(N+1,length(xp));
    
    % Initial values P_0(x) and P_1(x)
    if N==0
        %first component
        gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
        gamma(beta+1)/gamma(alpha+beta+1);
        PL(1,:) = 1.0/sqrt(gamma0);
        
    elseif N==1
        %first components
        gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
        gamma(beta+1)/gamma(alpha+beta+1);
        PL(1,:) = 1.0/sqrt(gamma0);
        
        %second components
        gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
        PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);

    else
        %first components
        gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
        gamma(beta+1)/gamma(alpha+beta+1);
        PL(1,:) = 1.0/sqrt(gamma0);
        
        %second components
        gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
        PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
        % Repeat value in recurrence.
        aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

        %other components
        % Forward recurrence using the symmetry of the recurrence.
        for i=1:N-1
            h1 = 2*i+alpha+beta;
            anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
            (i+1+beta)/(h1+1)/(h1+3));
            bnew = - (alpha^2-beta^2)/h1/(h1+2);
            PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
            aold =anew;
        end
    end
    PL = PL';
return