function [] = SetGlobalVariables(problemType, ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method,...
    Sylvester_Build_Method)
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

% Noise SEED for random numbers
SETTINGS.SEED = 1024;

% Outputs
SETTINGS.PLOT_GRAPHS = true;

%
% 'y' : Use Logs
% 'n'
%
SETTINGS.BOOL_LOG = false;

%--------------------------------------------------------------------------
%               
%                            Preprocessing
%
%
SETTINGS.BOOL_ALPHA_THETA = bool_alpha_theta;
SETTINGS.MEAN_METHOD = mean_method;

% -------------------------------------------------------------------------
% 
%                   SETTINGS IN COMPUTING GCD DEGREE
%

% Set the metric for measuring the degree of the GCD.
%
%   * Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals
%
SETTINGS.RANK_REVEALING_METRIC = 'Singular Values';


% Set the threshold for measuring the degree of the GCD. If max change in
% metric is less than this value, then all subresultants are full rank or
% rank deficient.
SETTINGS.THRESHOLD = 1;
SETTINGS.THRESHOLD_RANK = -5;

% -------------------------------------------------------------------------
%
%               SETTINGS FOR COMPUTING GCD COEFFICIENTS
%
%

% When computing the GCD do we wish to use both u(x,y) and v(x,y) or only
% u(x,y). Note: It is necessary to use both for standard GCD computations.
%
% 'ux and vx'
% 'ux'
%
SETTINGS.GCD_COEFFICIENT_METHOD = 'ux and vx';


%--------------------------------------------------------------------------

% Structuring the Sylvester Matrix
SETTINGS.SYLVESTER_BUILD_METHOD = Sylvester_Build_Method;

% T : 
% DT :
% DTQ :
% TQ : 
% DTQ Rearranged Denom Removed :
% DTQ Rearranged :




%-------------------------------------------------------------------------
%
%               LOW RANK SYLVESTER MATRIX SETTINGS.
%


% Regarding the computation of the Sylvester low rank approximation


% LOW_RANK_APPROXIMATION_METHOD:
%
% None : No perturbations added
% Standard STLN : Use linear form
% Root Specific SNTLN : Non-Linear low form, with g = f' constraint
% Standard SNTLN : Non-Linear low rank form, without constraint
%

SETTINGS.LOW_RANK_APPROXIMATION_METHOD = low_rank_approx_method;
SETTINGS.MAX_ERROR_SNTLN = 1e-13;
SETTINGS.MAX_ITERATIONS_SNTLN = 10;

%--------------------------------------------------------------------------
%
%           APF RELATED SETTINGS.
%
%

% Structuring the matrix [C(f) | C(g)]
SETTINGS.APF_METHOD = apf_method;

% APF_BUILD_METHOD
% * Standard
% * Rearranged
SETTINGS.APF_BUILD_METHOD = 'Standard';

% Regarding the computation of the low rank approximation of the 
% C = [C(u) ; C(v)] matrix
SETTINGS.MAX_ERROR_APF = 1e-14;
SETTINGS.MAX_ITERATIONS_APF = 50;

%-------------------------------------------------------------------------
%
%           INPUT VALIDATION
%
%


% Validation
% If BOOL_Q has not been included, then the Sylvester rearrangement is not
% applicable, and the common denominators can not be removed.
% Simplest method, no structure added.
% Override users input options if incompatable.


if ( strcmp(SETTINGS.PROBLEM_TYPE,'GCD') && strcmp(SETTINGS.LOW_RANK_APPROXIMATION_METHOD, 'Root Specific SNTLN'))
    
    fprintf([mfilename ' : Can not use root specific SNTLN method for GCD type problem']);
    SETTINGS.LOW_RANK_APPROXIMATION_METHOD = 'Standard SNTLN';
    
end
    

%



%-------------------------------------------------------------------------
%
%                   DECONVOLUTION SETTINGS
%
%
% Deconvolution method in the root finding problem 

%   'From Deconvolution'
%   'From ux'
SETTINGS.GET_HX_METHOD = 'From Deconvolution';

%
% Separate
% Batch
% Batch With STLN
% Batch Constrained
% Batch Constrained With STLN

deconvolve_method_hx = 'Batch';
SETTINGS.DECONVOLVE_METHOD_HX_FX = deconvolve_method_hx;

% Separate
% Batch
deconvolve_method_wx = 'Separate';
SETTINGS.DECONVOLVE_METHOD_WX_HX = deconvolve_method_wx;


SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-15;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 100;
SETTINGS.PREPROC_DECONVOLUTIONS = 'y';

end




