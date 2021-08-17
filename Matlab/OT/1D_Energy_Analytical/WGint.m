function [w] = WGint(alpbet,N,x)

    %Computes de gauss GLL quadrature weights for jacobi polynomials where
    % alpha = beta through the linear system of equations method
    % alpbet = alpha = beta
    % N = Polynomial jacobi family order
    % x = Jacobi polynomial (alpha = beta = alpbet) roots

    %Gauss weights as integrals
    K = ones(N+1,N+1);
    b = ones(N+1,1);
    for j=0:N
        K(j+1,:) = x.^j;
        %b(j+1) = ((1)^(j+1) - (-1)^(j+1))/(j+1);   %Legendre case
        %b(j+1) = ((-1)^j + 1)*gamma(alpbet+1)*gamma((j+1)/2)/(2*gamma((j+3)/2 + alpbet));
        b(j+1) = (((-1)^j + 1)*gamma(alpbet+1)*gamma((j+1)/2))/(2*gamma(alpbet + j/2 + 3/2));
    end
    w = K\b;
end