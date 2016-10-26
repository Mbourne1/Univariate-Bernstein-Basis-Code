function [  ] = Bisection_roots( example_number )
%ROOTS_BISECTION Summary of this function goes here
%   Detailed explanation goes here


addpath 'Examples'

[f_roots_exact,~] = Root_Examples(example_number);
f_exact_bi = BuildPolyFromRoots(f_roots_exact);

% Get the degree of f(x)
m = GetDegree(f_exact_bi);

% Get f(x) without binomial coefficients, in standard bernstein form.
f_exact = GetWithoutBinomials(f_exact_bi);

fx = f_exact;

% Get the set of control points of the curve
a = 0;
b = 1;

Pk = [];
for i = 0:1:m
    Pk = [Pk ; a + (i/m).*(b-a)     fx(i+1)];
end

% Let V be the set of control points.
Roots = []
FindRoots(Roots,Pk)



end

function [Roots] = FindRoots(Roots,V)


if (control_polygon_intesections == 0)
    fprintf('No roots');
elseif (control_polygon_intersections == 1)
    fprintf('control polygon is flat enough')
    fprintf('return t-intercept of chord control point 1 to n as the root')
    roots = [1 1];
else
    fprintf('subdivide curve at midpoint')
    % Subdivide at midpoint
    c = 0.5;
    [controlpoints_left,controlpoints_right] = deCasteljau(c,V)
    Roots = [Roots ; FindRoots(controlpoints_left)]
    Roots = [Roots ; FindRoots(controlpoints_right)]
end

end



function [new_left_pk, new_right_pk] = deCasteljau(c, Pk_0)
% Obtain the next set of control points, given a set of control points Pk_0
% in the deCasteljau subdivision process
% a = start point on horizontal axis
% b = end point on horizontal axis
% c = subdivision point on horizontal axis
% Pk_0 = original control points

% get first control poitn
a = Pk_0(1,1);
% get last control point
b = Pk_0(end,1);


T = (c-a)./(b-a);


% Build set of control points
Pk_1 = [];
% for
for i = 1:1:length(Pk_0)-1
    %Pk_1(i)
    aa =  (1-T).*Pk_0(i,:) + T.*Pk_0(i+1,:) ;
    Pk_1 = [Pk_1 ; aa];
    
end

% From the generated sets of control points obtain set of control points
% for sub-curve P[a,c]
new_left_pk = zeros(m+1,2);
for i = 1:1:m+1
    new_left_pk(i,:) = Pk_array{i}(1,:);
end

% obtain set of control points for sub-curve P[c,b]
new_right_pk = zeros(m+1,2);
for i=1:1:m+1
    new_right_pk(i,:) = Pk_array{i}(end,:);
end

end