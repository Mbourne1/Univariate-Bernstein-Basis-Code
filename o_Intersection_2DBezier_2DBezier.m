function [] =  o_Intersection_2DBezier_2DBezier()
% Given an example number, get the control points of 2 Bezier Curves, and
% calculate their intersections.
%
% 1. Get Control Points of polynomial f.
% 2. Get Control Points of polynomial g.
% 3. Implicitize Polynomial g.
% 4. Substitute the parametric expressions of f for x and y in t, into the
%    implicit polynomial g.

addpath('Root Finding Methods')


global PLOT_GRAPHS
global BOOL_Q
global BOOL_LOG
global BOOL_DENOM_SYL
global BOOL_PREPROC
global GEOMETRIC_MEAN_METHOD
global BOOL_SYLVESTER_BUILD_METHOD
global BOOL_SNTLN
global BOOL_APF
global BOOL_APF_BUILD_METHOD
global BOOL_DECONVOLVE
global NOMINAL_VALUE
global MIN_DELTA_MAG_ROW_SUMS
global BOOL_SNTLN_ROOTS
global MAX_ERROR
global MAX_ITERATIONS

PLOT_GRAPHS = 'y';
BOOL_Q = 'y';
BOOL_LOG = 'y';
BOOL_DENOM_SYL = 'y';
BOOL_PREPROC = 'y';
GEOMETRIC_MEAN_METHOD = 'matlab';
BOOL_SYLVESTER_BUILD_METHOD = 'standard';
BOOL_SNTLN = 'n';
BOOL_APF = 'n';
BOOL_APF_BUILD_METHOD = 'standard' ;
BOOL_DECONVOLVE = 'single';
BOOL_SNTLN_ROOTS = 'StandardSNTLN';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NOMINAL_VALUE = 10;
MAX_ERROR = 1e-15;
MAX_ITERATIONS = 40;

MIN_DELTA_MAG_ROW_SUMS = 2.5;

ex_num = '1';
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

Plot2DBezier(array_cp);


%% Implicitize the Bezier Control Points of g

fprintf('Implicit representation of g:\n')
[gxy,symbolic_expression] = ImplicitizeBezierBySylvester(CP_g);

PrintCoefficientsBivariate(gxy,'C_{2}:')

%%
% Get the parametric expressions of polynomial f in the brn basis
f_x =  CP_f(1,:)';
f_y =  CP_f(2,:);


% Put the sets of control points into an array
CP_arr{1} = CP_f;
CP_arr{2} = CP_g;


%% Substitute x(t) and y(t) into g

% get degree of input polynomial 
[~,c] = size(CP_f);
n = c-1;

[r,c] = size(gxy);

coef_Bernstein_Poly = zeros((n^2)+1,1);

% for each row in the implicit representation of g(x,y)
for i = 0:1:r-1
    % for each column
    for j = 0:1:c-1
        
        
        % Get the coefficient in g(x,y)
        coef = gxy(i+1,j+1);
        
        % Multiply x(t) from curve f by itself i times
        x_component = 1;
        for k = 1:1:i   
            x_component = BernsteinMultiply(x_component,f_x);
        end
        
        % Multiply y(t) from curve f by itself j times
        y_component = 1;
        for k = 1:1:j
            y_component = BernsteinMultiply(y_component,f_y');
        end
        
        % Muliply x(t)^i and y(t)^j
        xy_comp = BernsteinMultiply(x_component,y_component);
        
        
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

%% Get the roots of the polynomial


roots = o_roots_mymethod(coef_Bernstein_Poly);
roots = o_roots_matlab(coef_Bernstein_Poly)

%% 
% Given roots are calculated, plug in roots to one of the original 
% parametric equations
[num_rts,~] = size(roots);

% for each root
for i = 1:1:num_rts
    r_i = roots(i);
    if isreal(r_i)
    x = BernsteinEvaluate(f_x,r_i);
    y = BernsteinEvaluate(f_y',r_i);
    fprintf('The root %2.4e gives intersection point %2.3f , %2.3f \n',r_i,x,y)
   
    end
end
end

