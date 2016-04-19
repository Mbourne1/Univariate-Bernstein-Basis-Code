function [ root_mult_array ] = o_roots_subdivision(fx)
% Given the coefficients of the polynomial in Bernstein form, compute the
% roots using the subdivision method by Schneider in Graphics Gems.
%
% Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% Outputs.
%
% t : Roots

% Global variables
global MAXDEPTH
MAXDEPTH = 10;

global EPSILON
EPSILON = 1e-12;

% Initialise iteration number 
ite_num = 1;

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get the set of control points of the curve
a = 0;
b = 1;

% Initialise a matrix of control points, where each row contains the (x,y)
% pair.
Pk = zeros(m+1,2);

% Get the control points
for i = 0:1:m
    % Get the x ordinate
    x_ord = a + (i/m).*(b-a);
    Pk(i+1,:) = [x_ord fx(i+1)];
end

% Plot the coefficients and the control points of the curve
figure_name = sprintf('%s : Bernstein Polynomial Curve %i',mfilename,ite_num);
Plot_fx(fx,a,b,figure_name);

% Set the degree of polynomial f
degree = m;

W_DEGREE = m;

% Set the intial depth to 1
ini_depth = 1;

% Set the number of found roots
root_count = 0;

% set the current reached depth to 1
reached_depth = 1;

% Obtain the roots using the bisection method.
[t,~,reached_depth] = FindRoots(Pk,degree,ini_depth,W_DEGREE,root_count,reached_depth);

root_mult_array = [];

if isempty(t)
    fprintf('\nROOTS CALCULATED BY SUBDIVISION FUNCTION \n');
    fprintf('No Roots were Found\n')
    return
end

% %
% Given the set of calculated simple roots produce a polynomial in the 
% bernstein form.
roots = [real(t(:,1)), ones(length(t),1)];

% Get coefficients of polynomial r1 in the scaled bernstein form
r1_bi = B_poly(roots);

% Divide by binomial coefficients to obtain r1 in standard bernstein basis.
r1_x = GetWithoutBinomials(r1_bi);

% Obtain the polynomial f_{2}, given by the removal of the calculated roots
% from f_{1}.
f2 = Bernstein_Deconvolve(fx,r1_x);

% While f_{2} is not a scalar.
while GetDegree(f2) >= 1
    
    % Increment iteration number
    ite_num = ite_num + 1;
    
    % Get degree of polynomial f2
    m_ite = GetDegree(f2);
    
    % Get the control points of f2
    Pk = GetControlPoints(a,b,f2);
    
    % Plot f2
    figure_name = sprintf('%s : Bernstein Polynomial Curve %i',mfilename,ite_num);
    Plot_fx(f2,a,b,figure_name);
    
    % Get roots of polynomial f2
    [t_new,~,reached_depth] = FindRoots(Pk,m_ite,1,m_ite,1,1);
    
    % add the found roots to a list of roots
    root_mult_array = [root_mult_array; t_new ones(length(t_new))];
    
    % if no new roots are found, end this while loop
    [nEntries_t,~] = size(t_new);
    if (nEntries_t == 0)
        break;
    end
    
    % Build Polynomial from newly obtained roots
    roots = [real(t_new(:,1)), ones(length(t_new),1)];

    % Get polynomial coefficients of r1 in scaled bernstien basis
    r1_bi = B_poly(roots);
    
    % Get r1 in standard bernstein basis.
    r1_x = GetWithoutBinomials(r1_bi);
    
    % Perform deconvolution and obtain remaining part of original input
    % polynomial. Update f2 to be the result of the deconvolution.
    f2 = Bernstein_Deconvolve(f2,r1_x);
    
    
end

PrintoutRoots('SUBDIVISION' , root_mult_array);

end

function [roots,root_count,reached_depth] = FindRoots(CP,degree,curr_depth,W_DEGREE,root_count,reached_depth)
%
% Inputs
%
% CP - the set of control points [x,y]
%
% degree - degree of the polynomial
%
% curr_depth - current depth of iteration
%
% W_Degree - Degree of control points
%
% root count - Number of roots obtained to this depth.
%
% reached depth - maximum reached depth
%
% bool_printMessageLog - Print Messages.
%
% t -
%


% If the current depth is greater than reached depth so far, increment by 1
if curr_depth > reached_depth
    reached_depth = curr_depth;
end

global MAXDEPTH

% Initialise empty set of roots
roots = [];

% Get number of crossings of control polygon and x axis.
nCrossings = CrossingCount(CP);

% get the number of crossings of the control polygon.
switch(nCrossings)
    case 0 % No solutions in the interval
        
        
    case 1
        % If number of crossings is 1, then there is a unique solution in
        % the interval
        if (curr_depth >= MAXDEPTH) % Current depth is equal to or greater than the maximum depth
            
            % Calculate an approximate intercept, in the middle of the
            % interval.
            t_new = (CP(1,1) + CP(W_DEGREE+1,1)) ./2;
            
            % Add the approximate root to the set of roots
            roots = [roots;t_new];
            
            % Increment the root count.
            root_count = root_count + length(t_new);
            
           
            
        elseif (ControlPolygonFlatEnough(CP,degree))
            % if the control polygon is considered flat enough, that the
            % curve which it bounds may be considered a line. Calculate
            % the x intercept
            
            % Calculate the new root
            t_new = ComputeXIntercept(CP,degree);
            
            % Add the root to the list of roots.
            roots = [roots;t_new];
            
            % Increment the root count
            root_count = root_count + length(t_new);

            
        else
            % If the control polygon isn't flat enough. subdivide the
            % control polygon and look at the left and right control
            % polygons individually.
            
            % Get start point
            a = CP(1,1);
            
            % Get end point
            b = CP(degree+1,1);
            
            % Get midpoint
            c = a + ((b-a) /2);
            
            % Get a set of control points for the left and right halves of
            % the interval.
            [Left_Control_Points,Right_Control_Points] = BezierSubdivide(CP,degree,c);
            
            
            % Perform find roots on the left partition, while incrementing
            % the current depth by one.
            [t_left,~, reached_depth]  = FindRoots(Left_Control_Points,degree,curr_depth+1,W_DEGREE,root_count,reached_depth);
            
            % Add list of roots returned by findRoots() on the left
            % partition, to the list of all roots.
            roots = [roots;t_left];
            
            % Increment the root count.
            root_count = root_count + length(t_left);
            
            
            % Perform findRoots() on the right set of control points
            [t_right,~,reached_depth] = FindRoots(Right_Control_Points,degree,curr_depth+1,W_DEGREE,root_count,reached_depth);
            
            % Add list of roots returned by findRoots() on the right
            % partition, to the list of all roots.
            roots = [roots;t_right];
            
            % Increment root count.
            root_count = root_count + length(t_right);
            
        end
        
    otherwise
        % If the number of crossings is not 0 or 1, then num > 1.
        % The interval contains more than one root.
        
        if (curr_depth >= MAXDEPTH) % Current depth is equal to or greater than the maximum depth
            
            % Calculate an approximate intercept, in the middle of the
            % interval.
            t_new = (CP(1,1) + CP(W_DEGREE+1,1)) ./2;
            
            % Add the approximate root to the set of roots
            roots = [roots;t_new];
            
            % Increment the root count.
            root_count = root_count + length(t_new);
        else
        
        % Get Start Point
        a = CP(1,1);
        
        % Get End Point
        b = CP(degree+1,1);
        
        % Get Midpoint
        c = a + ((b-a) /2);
        
        % Get a set of control points for the left and right halves of
        % the interval.
        [Left_Control_Points,Right_Control_Points] = BezierSubdivide(CP,degree,c);
        

        % Perform find roots on the left partition, while incrementing
        % the current depth by one.
        [t_left,~,reached_depth]  = FindRoots(Left_Control_Points,degree,curr_depth+1,W_DEGREE,root_count,reached_depth);
        
        % Add list of roots returned by findRoots() on the left
        % partition, to the list of all roots.
        roots = [roots;t_left];
        
        % Increment root count.
        root_count = root_count + length(t_left);

        % Perform findRoots() on the right set of control points
        [t_right,~,reached_depth] = FindRoots(Right_Control_Points,degree,curr_depth+1,W_DEGREE,root_count,reached_depth);
        
        % Add any calculated roots to the list of all roots.
        roots = [roots;t_right];
        
        % Increment root count.
        root_count = root_count + length(t_right);
        end
end


end

function intercept = ComputeXIntercept(CP,degree)
% Calculate the intercept of the X axis of the line between the first and
% last control point.
%
% Inputs.
%
% CP :  Set of control points [x,y]
%
% degree :  Degree of control points.
%
%
%
%                       Outputs
%
% intercept :   Intercept of the x axis between C(0,0) and C(n,n)
%
%

% set first control point to be (x0,y0)
x_0 = CP(1,1);
y_0 = CP(1,2);

% set last control point to be (x1,y1)
x_1 = CP(degree+1,1);
y_1 = CP(degree+1,2);

% obtain change in y
d_y = y_1 - y_0;

% obtain change in x
d_x = x_1 - x_0;

% obtain gradient
m = d_y ./ d_x;

% calculate intercept of x axis by the line between (x0,y0) and (x1,y1)
intercept = x_1 - (y_1./m) ;

end

function n_crossings = CrossingCount(CP)
% Get number of crossings between a given control polygon and the X axis.
%
%                           Inputs
%
% CP :  Control Points [x,y]
%
% Outputs.
%
% n_crossings : Number of crossings
%


% Initialise old sign to be the value of y for the first control point.
old_sign = CP(1,2);

% Initialise number of crossings = 0.
n_crossings = 0;

% For each control point get its sign, and the sign of the previous control
% point, and compare them.
for i = 1:1:size(CP,1)-1
    
    % Assign the new sign to be y value for current control point
    new_sign = CP(i+1,2);
    
    % if there is a change of sign y0*y1 < 0, so incrememnt number of
    % crossings.
    if (new_sign * old_sign) < 0
        
        % if change of sign, then increment number of crossings
        n_crossings = n_crossings + 1;
        
    end
    
    % set old sign to be current iterations new sign, then start next
    % iteration.
    old_sign = new_sign;
end

end


function [new_left_pk, new_right_pk] = BezierSubdivide(Pk,degree,c)
% Given a set of control points, subdivide such that two new sets of
% control points are obtained.
%
% Inputs.
%
%
% Pk :  Set of original control points.
%
% degree :  degree of control points
%
% c :    point of subdivision.
%
%
% Outputs.
%
%
% new_left_pk : - The set of control points on the left of c.
%
% new_right_pk : - The set of control points on the right of c.
%


% obtain each set of control points
Pk_array{1} = Pk;

for i = 2:1:degree+1
    Pk_array{i} = deCasteljau(c,Pk_array{i-1})  ;
end

% From the generated sets of control points obtain set of control points
% for sub-curve P[a,c]
new_left_pk = zeros(degree+1,2);
for i = 1:1:degree+1
    new_left_pk(i,:) = Pk_array{i}(1,:);
end

% obtain set of control points for sub-curve P[c,b]
new_right_pk = zeros(degree+1,2);
for i=1:1:degree+1
    new_right_pk(i,:) = Pk_array{i}(end,:);
end

new_right_pk(i,:) = Pk_array{i}(end,:);

new_right_pk = sortrows(new_right_pk,1);

end

function [Pk_1] = deCasteljau(c, Pk_0)
% Obtain the next set of control points, given a set of control points Pk_0
% in the deCasteljau subdivision process.
% Uses deCasteljau algorithm found in 'Applied Geometry for computer
% graphics and CAD' - Page 151

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



end


function val = ControlPolygonFlatEnough(CP,degree)
% Given a control polygon, return a boolean as to whether it is flat enough
% to be considered a straight line.
%
% Inputs
%
%
% CP :  Set of control points which form the control polygon [x,y]
%
% degree :  Degree of control polygon
%
% Outputs
%
%
% val - (Boolean)
%   1 : Control polygon is flat enough to be considered a line
%   0 : Control polygon is not flat enough.


global EPSILON


% derive implicit equation fo line connect first and last point
a = CP(1,1) - CP(degree,1);
b = CP(1,2) - CP(degree,1);
c = (CP(1,1) * CP(degree,2)) - (CP(degree,1) * CP(1,2));

abSquared = (a.*a) + (b.*b);

for i = 1:1:degree
    distance(i) = a.*CP(i,1) + b.*CP(i,2) + c;
    if (distance > 0)
        distance(i) = (distance(i).*distance(i)) ./ abSquared;
    elseif (distance < 0)
        distance(i) = -((distance(i).* distance(i)))./ abSquared;
    end
end

max_distance_above = 0;
max_distance_below = 0;

for i =1:1:degree-1
    if (distance(i) < 0)
        fprintf('')
        max_distance_below = min(max_distance_below,distance(i));
    elseif (distance(i) > 0)
        max_distance_above = max(max_distance_above,distance(i));
    end
end

% Implicit equation for zero line
a1 = 0.0;
b1 = 1.0;
c1 = 0.0;

% implicit equation for "above" line
a2 = a;
b2 = b;
c2 = c+ max_distance_above;

det = (a1 .* b2) - (a2 .* b1);
dInv = 1.0 ./ det;

intercept_1 = (b1 * c2 - b2 * c1) * dInv;

% Implicit equation for "below" line

a2 = a;
b2 = b;
c2 = c+ max_distance_below;

det = (a1 .* b2) - (a2 .* b1);
dInv = 1.0 ./ det;

intercept_2 = (b1 * c2 - b2 * c1) * dInv;


% compute Intercepts of bounding box
left_intercept = min(intercept_1,intercept_2);
right_intercept = max(intercept_1,intercept_2);

error = 0.5 * (right_intercept-left_intercept);
if (error < EPSILON)
    %fprintf('Control Polygon is flat enough. \n')
    val = 1;
    
else
    %fprintf('Control Polygon is not flat enough. \n')
    val = 0;
end
end





