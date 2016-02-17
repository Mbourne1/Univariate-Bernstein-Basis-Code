function [] = o_intersection_Bernstein_Bernstein(ex_num)
% Given two Bernstein polynomials, calculate the points of intersection

[C1,C2] = Examples_Intersection(ex_num)

Bezier{1} = C1;
Bezier{2} = C2;

% Plot the curves
Plot_Bezier(Bezier);

% get f(x) and g(x)
fx = C1(2,:);
gx = C2(2,:);

% Add Noise

hx = fx - gx;

%% Get Roots by Matlab Method
o_roots_matlab(hx')

%% Get Roots by my Method
% Set global parameters for mymethod
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
bool_sntln = 'y';
bool_apf = 'n';
Bool_APFBuildMethod = 'standard' ;
bool_deconvolve = 'single';
bool_SNTLN_Roots = 'StandardSNTLN';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nominal_value = 5;
max_error = 1e-15;
max_iterations = 40;

min_delta_mag_rowsum = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o_roots_mymethod(hx')
%%
% Given the points of intersection get the points of intersection.
% Evaluate the Bézier curve at the points of intersection to find y
% coordinate.

end
