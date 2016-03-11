function [ft] = BernsteinEvaluate(f,t)
% Given the function f(t), evaluate f(t) a point t.
%
%   Input.
%   
%   f   :   Coefficients of polynomial f(t) in Bernstein form.
%
%   t   :   Point at which to evaluate f(t)

[r,~] = size(f);
m = r-1;

sum = 0;

for i = 0:1:m
    
    sum = sum + f(i+1) .* nchoosek(m,i) .* ((1-t).^(m-i)) .* (t^i);
    
end

ft = sum;

end


