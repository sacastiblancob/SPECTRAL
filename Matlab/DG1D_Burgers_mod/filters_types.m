%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proofs over filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vandeven filter for minimizing A
p = 50;
x = 0:0.01:1;
y = ones(length(x),1);
for i=1:length(x)
  a = x(i);
  y(i) = -((a^(-p))*((-(a - 1)*a)^p)*hypergeom([1-p p],p+1,1-a));
  y(i) = y(i)/p;
end
y = 0 - (factorial(2*p - 1)/(factorial(p - 1)^2))*y;
y(1) = 1.0;
plot(x,y)






