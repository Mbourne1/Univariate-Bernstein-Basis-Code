function Pk = GetControlPoints(a,b,f)
%% Get set of control points for the function f on the interval [a,b]
% a :- Interval Lower Limit.
% b :- Interval Upper Limit.
% f :- Coefficients of Polynomial
%
% Pk :- Set of control points.

% Get degree of polynomial f
m = length(f)-1;

% Initialise empty vector of control points.
Pk = [];

% for each control point, assign value.
for i = 0:1:m
    Pk = [Pk; a+(i/m).*(b-a)    f(i+1)];
end



end