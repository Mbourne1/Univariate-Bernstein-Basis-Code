function [] = SetGlobalVariables_Roots(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method,...
    Sylvester_Build_Method, rank_revealing_metric, deconvolution_method_hx, deconvolution_method_wx)
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
%
%
%
%
%
%

global SETTINGS

% Set example number
SETTINGS.EX_NUM = ex_num;

% Set noise levels
SETTINGS.EMIN = emin;
SETTINGS.EMAX = emax;


% Noise SEED for random numbers
SETTINGS.SEED = 1024;

% Outputs
SETTINGS.PLOT_GRAPHS = false;

% BOOL_LOG (Boolean)
%   true : Use Logs
%   false : Dont use logs
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
%   * Max/Min Singular Values
%   * R1 Row Norms
%   * R1 Row Diagonals
%   * Residuals
%
SETTINGS.RANK_REVEALING_METRIC = rank_revealing_metric;


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
% DTQ Denominator Removed :
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

% APF_BUILD_METHOD (String) 
% * Standard
% * Rearranged
SETTINGS.APF_BUILD_METHOD = 'Standard';

% Regarding the computation of the low rank approximation of the 
% C = [C(u) ; C(v)] matrix
SETTINGS.MAX_ERROR_APF = 1e-14;
SETTINGS.MAX_ITERATIONS_APF = 50;



%-------------------------------------------------------------------------
%
%                   DECONVOLUTION SETTINGS
%
%
% Deconvolution method in the root finding problem 

% GET_HX_METHOD
%   'From Deconvolution'
%   'From ux'
SETTINGS.GET_HX_METHOD = 'From Deconvolution';

% DECONVOLUTION_METHOD_HX_FX
% * Separate
% * Batch
% * Batch With STLN
% * Batch Constrained
% * Batch Constrained With STLN
SETTINGS.DECONVOLUTION_METHOD_HX = deconvolution_method_hx;

% DECONVOLUTION_METHOD_WX_HX
% Separate
% Batch
if (strcmp(deconvolution_method_wx, 'Batch') || strcmp(deconvolution_method_wx ,'Separate'))
    SETTINGS.DECONVOLUTION_METHOD_WX = deconvolution_method_wx;
else
    error('not valid')
end

SETTINGS.MAX_ERROR_DECONVOLUTIONS = 1e-15;
SETTINGS.MAX_ITERATIONS_DECONVOLUTIONS = 10;
SETTINGS.PREPROC_DECONVOLUTIONS = true;

end




