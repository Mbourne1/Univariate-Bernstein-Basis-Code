function gx = Bernstein_Differentiate(fx)
% Differentiate the polynomial f(x) whose coefficients are given in the
% Bernstein basis. Let the derivative of f(x) be known as g(x).
%
% Inputs.
%
% fx : Coefficients of input polynomial f(x)
%
% Outputs.
%
% gx : Deriveative of polynomial f(x).

% Get degree of polynomial f(x)
m = size(fx,1) - 1;

% Initialise the vector of the derivative of f
gx = zeros(m,1);

%Loop through all coefficients of g(x)
for i = 1:1:m
    gx(i+1) = m * (fx(i+1)-fx(i));
end




end