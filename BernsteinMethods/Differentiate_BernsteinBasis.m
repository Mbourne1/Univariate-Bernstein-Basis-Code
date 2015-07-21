function der_f = Differentiate_BernsteinBasis(f)
%% Differentiate the polynomial coefficients of f, in the Bernstein basis.
% f = coefficients of input polynomial.
% der_f = deriveative of polynomial f.

%Get degree of polynomial 
    m = length(f)-1;
    
% Initialise the vector of the derivative of f    
    der_f = zeros(1,m);
    
% Loop through all coefficients     
    for i = 0:1:m-1
        k = i+1;
        der_f(k) = m* (f(k+1)-f(k));
    end


end