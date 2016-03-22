function [px] = Bernstein_Multiply(fx,gx)
% Given the Coefficients of two polynomials f(x) and g(x) in vector form, 
% output the coefficients of the product p(x).
%
% % Input.
%
% fx : Column vector of coefficients of the Bernstein Polynomial f(x)
%
% gx : Column vector of coefficients of the Bernstein Polynomial g(x)
%
% % Output.
%
% px : Column vector of coefficients of the Bernstein Polynomial p(x)


% Get the degree of polynomial f(x)
m = size(fx,1) - 1;

% Get the degree of polynomial g(x)
n = size(gx,1) - 1;

% Binomial coefficients corresponding to f(x)
Bi_m = GetBinomials(m);

% Binomail coefficients corresponding to g(x)
Bi_n = GetBinomials(n);

% Get fw
f_bi = fx .* Bi_m;


% Build matrix C
% For each column, insert the coefficients of f(x) and corresponding
% binomial coefficient a_{i}\binom{m}{i}.
C = zeros(n+1,m+n+1);
for i = 0:1:n
    C(i+1:(m+1)+i,i+1) = f_bi;
end

% Build diagonal matrix D^{-1}
D = diag(1./GetBinomials(m+n));


% Build matrix Q
Q = diag(Bi_n);

% Get the vector of coefficients of p(x) from the product f(x) and g(x).
px = D*C*Q*gx;

end