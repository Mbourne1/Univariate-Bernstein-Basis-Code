function [] =  o_Intersection_2DBezier_2DBezier(ex_num,bool_preproc,low_rank_approx_method,apf_method)
% o_Intersection_2DBezier_2DBezier
%       ...(ex_num,bool_preproc,low_rank_approx_methed,apf_method)
%
% Given an example number, get the control points of 2 Bezier Curves, and
% calculate their intersections.
%
% Inputs.
% 
% ex_num :
%
% bool_preproc :
%
% low_rank_approx_method : 
%
% apf_method :
%
%


% 1. Get Control Points of polynomial f.
% 2. Get Control Points of polynomial g.
% 3. Implicitize Polynomial g.
% 4. Substitute the parametric expressions of f for x and y in t, into the
%    implicit polynomial g.

addpath('Root Finding Methods')


SetGlobalVariables(bool_preproc,low_rank_approx_method,apf_method);

switch ex_num
    case '1'
        ex_num_f = '1';
        ex_num_g = '1';
end


% Get a set of control points for the planar Bezier curve f
CP_f = Examples_2DBezier_ControlPoints(ex_num_f);

% Get a set of control points for the planar Bezier curve g
CP_g = Examples_2DBezier_ControlPoints(ex_num_g);

array_cp{1} = CP_f;
array_cp{2} = CP_g;

Plot2DBeziers(array_cp);


%% Implicitize the Bezier Control Points of g

fprintf('Implicit representation of g:\n')
[gxy,symbolic_expression] = ImplicitizeBezierBySylvester(CP_g);

% Print the coefficients of the polynomial g(x,y) in Bernstein form.
PrintCoefficients_Bivariate_Bernstein(gxy,'C_{2}:')

%%
% Get the parametric expressions of polynomial f in the brn basis



% Put the sets of control points into an array
CP_arr{1} = CP_f;
CP_arr{2} = CP_g;

% Get coefficients of f(x)
f_x =  CP_f(1,:)';

% Get Coefficients of f(y)
f_y =  CP_f(2,:);

% Get Coefficients of the Bernstein polynomial obtained by susbstituting
% x(t) and y(t) from f(x,y) into g(x,y).
coef_Bernstein_Poly = Substitute(CP_f,gxy);

%% Get the roots of the polynomial
roots = o_roots_mymethod(coef_Bernstein_Poly);
roots = o_roots_matlab(coef_Bernstein_Poly);

%% 
% Given roots are calculated, plug in roots to one of the original 
% parametric equations
[num_rts,~] = size(roots);

% for each root
for i = 1:1:num_rts
    r_i = roots(i);
    if isreal(r_i)
    x = Bernstein_Evaluate(f_x,r_i);
    y = Bernstein_Evaluate(f_y',r_i);
    fprintf('The root %2.4e gives intersection point %2.3f , %2.3f \n',r_i,x,y)
   
    end
end
end


function [coef_Bernstein_Poly]= Substitute(CP_f,gxy)
%% Substitute x(t) and y(t) from f(x,y) into g(x,y)

% Get degree of input polynomial 
[~,c] = size(CP_f);
n = c-1;

% Get coefficients of f(x)
f_x =  CP_f(1,:)';

% Get Coefficients of f(y)
f_y =  CP_f(2,:);

% Get Degree of g(x,y)
[n1,n2] = GetDegree_Bivariate(gxy);

% Initialise the univariate Bernstein basis polynomial obtained from the
% substitution.
coef_Bernstein_Poly = zeros((n^2)+1,1);

% for each row in the implicit representation of g(x,y)
for i = 0:1:n1
    % for each column of g(x,y)
    for j = 0:1:n2
        
        
        % Get the coefficient b_{i,j} in g(x,y)
        coef = gxy(i+1,j+1);
        
        % Coefficient b_{i,j} has x^{i}
        
        % Get x(t)^{i}
        x_component = 1;
        for k = 1:1:i   
            x_component = Bernstein_Multiply(x_component,f_x);
        end
        
        % Get y(t)^{j}
        y_component = 1;
        for k = 1:1:j
            y_component = Bernstein_Multiply(y_component',f_y')';
        end
        
        % Get x(t)^i * y(t)^j
        xy_comp = Bernstein_Multiply(x_component,y_component);
        
        
        uij =  xy_comp;
        
        % output polynomial will be of degree 2n
        m = n^2;
        
        % degree elevate uij
        [r1,~] = size(uij);
        curr_deg_uij = r1-1;
        num_deg_elv_req =  m-curr_deg_uij;

        uij = DegreeElevate_Univariate(uij,num_deg_elv_req);
        
        uij = coef .* uij;
        
        coef_Bernstein_Poly = coef_Bernstein_Poly + uij;
    end
end
end
