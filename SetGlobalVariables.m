function [] = SetGlobalVariables(bool_preproc,low_rank_approx_method, apf_method)
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


%% Structuring the Sylvester Matrix

global BOOL_Q
BOOL_Q = 'y';

global BOOL_DENOM_SYL
BOOL_DENOM_SYL = 'y';

global SYLVESTER_BUILD_METHOD
% Standard D*T*Q
% Rearranged : Remove common denominator
SYLVESTER_BUILD_METHOD = 'Rearranged';


%% Structuring the matrix [C(f) | C(g)]

global APF_METHOD
APF_METHOD = apf_method;

global BOOL_DENOM_APF
BOOL_DENOM_APF = 'y';

global APF_BUILD_METHOD
APF_BUILD_METHOD = 'Standard';


%% Preprocessing

global BOOL_PREPROC
BOOL_PREPROC = bool_preproc;

global GEOMETRIC_MEAN_METHOD
GEOMETRIC_MEAN_METHOD = 'matlab';

%% Numerical Considerations

global BOOL_LOG
BOOL_LOG = 'y';

%% Noise 

global SEED
SEED = 1024;

%% Outputs

global PLOT_GRAPHS
PLOT_GRAPHS = 'y';

%% Regarding the computation of the Sylvester low rank approximation

% LOW_RANK_APPROXIMATION_METHOD:
% None : No perturbations added
% Standard STLN : Use linear form
% Root Specific SNTLN : Non-Linear low form, with g = f' constraint
% Standard SNTLN : Non-Linear low rank form, without constraint

global LOW_RANK_APPROXIMATION_METHOD
LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;

global MAX_ERROR_SNTLN
MAX_ERROR_SNTLN = 1e-15;

global MAX_ITERATIONS_SNTLN
MAX_ITERATIONS_SNTLN = 50;

%% Regarding the computation of the low rank approximation of the 
% C = [C(u) ; C(v)] matrix

global MAX_ERROR_APF
MAX_ERROR_APF = 1e-12;

global MAX_ITERATIONS_APF
MAX_ITERATIONS_APF = 50;

% global BOOL_APF_ROOTS
% BOOL_APF_ROOTS = 'StandardAPF';


% nominal value used when:
% Only one sylvester subresultant exists, ie k = 1 and min(m,n) = 1. where
% m and n are the degrees of input polynomials f and g.
% if max_r./min_r > nominal_value (then minimum value is significantly
% small, to assume that the sylvester matrix is rank deficient)
% then degree is one. otherwise degree is zero
global THRESHOLD
THRESHOLD = 2;


%% Validation
% If BOOL_Q has not been included, then the Sylvester rearrangement is not
% applicable, and the common denominators can not be removed.
% Simplest method, no structure added.
% Override users input options if incompatable.
if (BOOL_Q == 'n')
    LOW_RANK_APPROXIMATION_METHOD = 'None';
    APF_METHOD = 'None'; % Does not work with code block APF (Addition of structured perturbation code doesnt exist for exclusion of Q from coefficient matrix).
    BOOL_DENOM_SYL = 'y';
end

%%
% Intialise the global variables
global MAX_ERROR_DECONVOLUTIONS
global MAX_ITERATIONS_DECONVOLUTIONS
global BOOL_DECONVOLVE

% Get Global variables for Deconvolve
MAX_ERROR_DECONVOLUTIONS = 1e-10;
MAX_ITERATIONS_DECONVOLUTIONS = 50;
BOOL_DECONVOLVE = 'Single';

