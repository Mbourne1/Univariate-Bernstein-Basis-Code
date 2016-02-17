function [] = PrintPoly(fx,name)
% Given polynomial coefficients in bernstein form, print in latex format
%
%   Inputs
%
%   fx : Coefficients of polynomial f(x) in the Bernstein form.
%
%   name : name of the function eg 'f'
%
%


% Get the degree of polynomial f
[r,~] = size(fx);
m = r - 1;

str = sprintf(' %s(x) = ',name);

for i = 0:1:m
    % Get the coefficient
    coeff = fx(i+1);
    
    % Create string
    if coeff>=0
        temp_str = sprintf('+%2.4f B_{%i}^{%i}', coeff, i, m);
    else
        temp_str = sprintf('%2.4f B_{%i}^{%i}', coeff, i, m);
    end
    % Append to string
    str = strcat(str,temp_str);
    
end

str = strcat(str,'\n');

% Print out the polynomial
fprintf(str);

end