function [ft] = BernsteinEvaluate(f,t)
% Given the function f(t), evaluate at a point t.

[r,~] = size(f);
m = r-1;

sum = 0;

for i = 0:1:m
    
    sum = sum + f(i+1) .* nchoosek(m,i) .* ((1-t).^(m-i)) .* (t^i);
    
end

ft = sum;

end


