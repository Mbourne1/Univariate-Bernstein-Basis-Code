function [c] = DegreeElevate(fx,r)
% Function performs degree elevation on polynomial f(x).
% Let fx be the coefficients of polynomial f(x) of degree m. 
% r = Number of degree elevations such that the output polynomial is of
% degree m + r.

n = length(fx) -1;

% for each c(i) in the vector
for k = 0:1:n+r
   temp_sum = 0;
   for j = max(0,k-r):1:min(n,k)
       temp_sum = temp_sum + (fx(j+1) * nchoosek(n,j) * nchoosek(r,k-j)) ./ nchoosek(n+r,k) ;
   end
   c(k+1) = temp_sum;
end