function int_f= Integrate(f)
% Integrate a polynomial given in the Bernstein Basis

% f = coefficients of input polynomial

% int_f = coefficients of integral of f.

% Get Degree of input polynomial.
n = length(f)-1;

% Initialise vector of integrated polynomial coefficients.
int_f = zeros(1,n+1);

% Loop through each coefficient

for k = 0:1:n+1
    i = k+1;
    if k > 0
        I = (1/(n+1)) * ...
            sum(f(1:(i-1)));
    else
        I = 0;
    end
    int_f(i) = I;
end


end