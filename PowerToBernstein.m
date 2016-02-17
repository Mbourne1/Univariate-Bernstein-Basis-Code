function [C] = PowerToBernstein(fx)
% given a set of coefficients in the power basis, convert to bernstein
% basis
% Algorithms for Polynomials in Bernstein Form


% let n be the degree of polynomial fx
n = length(fx)-1;


for j = 0:1:n
    temp_sum = 0;
    for k = 0:1:j
        temp_sum = temp_sum + ...
            (...
                nchoosek(j,k) ./ nchoosek(n,k) .* fx(k+1)...
            );
    end
    C(j+1) = temp_sum;
end
