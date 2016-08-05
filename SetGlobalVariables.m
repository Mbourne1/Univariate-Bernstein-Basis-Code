function [] = SetGlobalVariables(problemType, ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
% Set the global variables
%
% Inputs.
%
% bool_preproc
%
% low_rank_approx_method
%
% apf_method
%

global SETTINGS
SETTINGS.EX_NUM = ex_num;

SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;

% Set the problem Type
SETTINGS.PROBLEM_TYPE = problemType;

%
% 'y' : include the matrix Q in the Sylvester matrix S = D^{-1}T(f,g)Q
% 'n' : Sylvester matrix excludes Q, S = D^{-1}T(f,g)
%
SETTINGS.BOOL_Q = 'y';

%
% 'y' : Use Logs
% 'n'
%
SETTINGS.BOOL_LOG = 'n';

% When computing the GCD do we wish to use both u(x,y) and v(x,y) or only
% u(x,y). Note: It is necessary to use both for standard GCD computations.
%
% 'ux and vx'
% 'ux'
%
SETTINGS.GCD_COEFFICIENT_METHOD = 'ux and vx';


% Get GCD Degree

% Set the metric for measuring the degree of the GCD.
SETTINGS.METRIC = 'Singular Values';

% Set the threshold for measuring the degree of the GCD. If max change in
% metric is less than this value, then all subresultants are full rank or
% rank deficient.
SETTINGS.THRESHOLD = 1;
SETTINGS.THRESHOLD_RANK = -5;

% Structuring the Sylvester Matrix
SETTINGS.BOOL_DENOM_SYL = 'y';
SETTINGS.SYLVESTER_BUILD_METHOD = 'Standard';

% Standard : D*T*Q
%
% Rearranged : In a format where the entries of the sylvester matrix
% partition C(f) have a common divisor, and the entries of C(g) have a
% common divisor


% Structuring the matrix [C(f) | C(g)]
SETTINGS.APF_METHOD = apf_method;
SETTINGS.BOOL_DENOM_APF = 'y';
SETTINGS.APF_BUILD_METHOD = 'Rearranged';


% Preprocessing

SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;
SETTINGS.MEAN_METHOD = mean_method;


% %
% %
% Numerical Considerations


% Noise 
SETTINGS.SEED = 1024;

% Outputs
SETTINGS.PLOT_GRAPHS = 'y';

% Regarding the computation of the Sylvester low rank approximation

% LOW_RANK_APPROXIMATION_METHOD:
% None : No perturbations added
% Standard STLN : Use linear form
% Root Specific SNTLN : Non-Linear low form, with g = f' constraint
% Standard SNTLN : Non-Linear low rank form, without constraint

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;
SETTINGS.MAX_ERROR_SNTLN = 1e-16;
SETTINGS.MAX_ITERATIONS_SNTLN = 75;

% Regarding the computation of the low rank approximation of the 
% C = [C(u) ; C(v)] matrix
SETTINGS.MAX_ERROR_APF = 1e-14;
SETTINGS.MAX_ITERATIONS_APF = 50;




% Validation
% If BOOL_Q has not been included, then the Sylvester rearrangement is not
% applicable, and the common denominators can not be removed.
% Simplest method, no structure added.
% Override users input options if incompatable.



if (SETTINGS.BOOL_Q == 'n')
    SETTINGS.SYLVESTER_BUILD_METHOD = 'Standard';
    SETTINGS.LOW_RANK_APPROXIMATION_METHOD = 'None';
    SETTINGS.APF_METHOD = 'None'; % Does not work with code block APF (Addition of structured perturbation code doesnt exist for exclusion of Q from coefficient matrix).
    SETTINGS.BOOL_DENOM_SYL = 'y';
end


if ( strcmp(SETTINGS.PROBLEM_TYPE,'GCD') && strcmp(SETTINGS.LOW_RANK_APPROXIMATION_METHOD, 'Root Specific SNTLN'))
    
    fprintf([mfilename ' : Can not use root specific SNTLN method for GCD type problem']);
    SETTINGS.LOW_RANK_APPROXIMATION_METHOD = 'Standard SNTLN';
    
end
    

%
% Get Global variables for Deconvolve
SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-12;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 50;



