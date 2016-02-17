function [] =  o_intersection_Bezier_Bezier()
% Given an example number, get the control points of 2 Bezier Curves, and
% calculate their intersections.
%
% 1. Get Control Points of polynomial f.
% 2. Get Control Points of polynomial g.
% 3. Implicitize Polynomial g.
% 4. Substitute the parametric expressions of f for x and y in t, into the
%    implicit polynomial g.

addpath('Root Finding Methods')
addpath('BernsteinMethods')

global bool_plotgraphs
global bool_q
global bool_log
global bool_denom_syl
global bool_preproc
global geometricMeanMethod
global bool_sylvesterBuildMethod
global bool_sntln
global bool_apf
global Bool_APFBuildMethod
global bool_deconvolve
global nominal_value
global min_delta_mag_rowsum
global bool_SNTLN_Roots
global max_error
global max_iterations

bool_plotgraphs = 'y';
bool_q = 'y';
bool_log = 'y';
bool_denom_syl = 'y';
bool_preproc = 'y';
geometricMeanMethod = 'matlab';
bool_sylvesterBuildMethod = 'standard';
bool_sntln = 'n';
bool_apf = 'n';
Bool_APFBuildMethod = 'standard' ;
bool_deconvolve = 'single';
bool_SNTLN_Roots = 'StandardSNTLN';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nominal_value = 10;
max_error = 1e-15;
max_iterations = 40;

min_delta_mag_rowsum = 2.5;

ex_num = '1'
switch ex_num
    case '1'
        ex_num_f = '1';
        ex_num_g = '2';
end


% Get a set of control points for a polynomial f
CP_f = Examples_Bezier_ControlPoints(ex_num_f);

% Get the equation of the second bezier curve
CP_g = Examples_Bezier_ControlPoints(ex_num_g);

array_cp{1} = CP_f
array_cp{2} = CP_g

Plot_Bezier(array_cp)


%% Implicitize the Bezier Control Points of g

fprintf('Implicit representation of g:\n')
[implicit_g,symbolic_expression] = Implicitize_Bezier_Sylvester(CP_g)


%%
% Get the parametric expressions of polynomial f in the brn basis
f_x =  CP_f(1,:)';
f_y =  CP_f(2,:);


% Put the sets of control points into an array
CP_arr{1} = CP_f;
CP_arr{2} = CP_g;


%% Substitute x(t) and y(t) into g

% get degree of input polynomial 
[~,c] = size(CP_f)
n = c-1;

[r,c] = size(implicit_g);

coef_Bernstein_Poly = zeros( (n^2)+1,1);

% for each row in the implicit representation of g(x,y)
for i = 0:1:r-1
    % for each column
    for j = 0:1:c-1
        
        
        % Get the coefficient in g(x,y)
        coef = implicit_g(i+1,j+1);
        
        % Multiply x(t) from curve f by itself i times
        x_component = 1;
        for k = 1:1:i   
            x_component = Bernstein_Multiply(x_component,f_x);
        end
        
        % Multiply y(t) from curve f by itself j times
        y_component = 1;
        for k = 1:1:j
            y_component = Bernstein_Multiply(y_component,f_y');
        end
        
        % Muliply x(t)^i and y(t)^j
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
        
        coef_Bernstein_Poly = coef_Bernstein_Poly + uij
    end
end

%% Get the roots of the polynomial


roots = o_roots_mymethod(coef_Bernstein_Poly)
roots2 = o_roots_matlab(coef_Bernstein_Poly)

%% Given roots are calculated, plug in roots to one of the original parametric equations
[num_rts,~] = size(roots)

% for each root
for i = 1:1:num_rts
    r_i = roots(i)
    if isreal(r_i)
    x = BernsteinEvaluate(f_x,r_i);
    y = BernsteinEvaluate(f_y',r_i);
    fprintf('Intersection Point: \n [%2.3f , %2.3f] \n',x,y)
    end
end
end

