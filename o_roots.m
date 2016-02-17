function [] = o_roots(ex_num,emin,emax,BOOL_SNTLN,BOOL_APF,BOOL_PREPROC,seed)
% Given an example number, and a set of input parameters, calculate the
% roots r_{i} of the polynomial f(x) and the corresponding multiplicities 
% m_{i} 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           Inputs
%
% ex - (Int) Example Number
%
% emin - Noise/Signal maximum threshold (minimum)
%
% emax - Noise/Signal maximum threshold (maximum)
%
% BOOL_SNTLN - Assigned to global variable (see below)
%
% BOOL_APF - Assigned to global variable (see below)
%
% BOOL_DENOM - Assigned to global variable (see below)
%
% BOOL_PREPROC - Assigned to global variable (see below)
%
% seed (int) - Integer chosen to randomly generate the roots and
% multiplicities of input polynomials f and g as well as the noise which is
% added to their coefficients.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                    # Global Variables #


global bool_bezout
global bool_sntln
global bool_apf
global bool_denom_syl
global bool_denom_apf
global bool_preproc
global bool_q
global bool_log

global PLOT_GRAPHS

global MAX_ERROR_SNTLN
global MAX_ITERATIONS_SNTLN

global MAX_ERROR_APF
global MAX_ITERATIONS_APF

global SEED

global bool_sylvesterBuildMethod
global Bool_APFBuildMethod
global max_error_deconvolutions
global max_iterations_deconvolutions
global problemType 
global bool_deconvolve

bool_sntln = BOOL_SNTLN;
bool_preproc = BOOL_PREPROC;
bool_apf = BOOL_APF;
bool_bezout = 1;
bool_denom_syl = 'y';
bool_denom_apf = 'y';
bool_q = 'y';
bool_log = 'y';
PLOT_GRAPHS = 'n';
SEED = seed;

MAX_ERROR_SNTLN = 1e-15;
MAX_ITERATIONS_SNTLN = 50;

MAX_ERROR_APF = 1e-15;
MAX_ITERATIONS_APF = 50;

max_error_deconvolutions = 1e-10;
max_iterations_deconvolutions = 50;

bool_sylvesterBuildMethod = 'rearranged';
Bool_APFBuildMethod = 'rearranged';

problemType = 'fromRoots'; % fromRoots/fromCoefficients
bool_deconvolve = 'batch';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nominal value used when:
% Only one sylvester subresultant exists, ie k = 1 and min(m,n) = 1. where
% m and n are the degrees of input polynomials f and g.
% if max_r./min_r > nominal_value (then minimum value is significantly
% small, to assume that the sylvester matrix is rank deficient)
% then degree is one. otherwise degree is zero
global nominal_value
nominal_value = 10;

% let x be the maximum change in ratio_maxmin_rowsum vector if abs(x) <
% nominal_value_2, if the change is minimal, then all subresultants should
% be classed as rank deficient.

global min_delta_mag_rowsum
min_delta_mag_rowsum = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bool_SNTLN_ROOTS
% RootSpecificSNTLN :   Use Roots based SNTLN method, which has the added constraints that
%                       g is the derivative of f.
% StandardSNTLN :   Use standard SNTLN where f and g are unconstrained
global bool_SNTLN_Roots
bool_SNTLN_Roots = 'StandardSNTLN';

% bool_APF_Roots
% RootSpecificAPF :   Use roots based APF method, which has added constraings.
% StandardAPF :   Use standard apf method where f and g are unconstrained.
global bool_APF_Roots
bool_APF_Roots = 'StandardAPF';


% geometricMeanMethod
% used when calculating geometric means of the entries of a Sylvester
% matrix, a standard method would consider each entry in the matrix, but a
% new method, described in my internal report offers speed up due to the
% structured nature of the Sylvester matrix entries for the Sylvester
% matrix in the Bernstein basis.
% 'matlab' -   use MatLab Built in method for calculating geometric means
% 'mymethod' -   use my method of calculating geometric means#
% 'none' -   Set geometric mean equal to one, essentially eliminating the
%       preprocessor.
global geometricMeanMethod
geometricMeanMethod = 'matlab';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Validate Inputs.


% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
if emin > emax
    fprintf('minimum noise greater than maximum noise \n swapping values...\n')
    emin_wrong = emin;
    emax_wrong = emax;
    emin = emax_wrong;
    emax = emin_wrong;
end

% If BOOL_Q has not been included, then the Sylvester rearrangement is not
% allowed, and the denominator can not be removed.
switch bool_q
    case 'n'
        
        bool_denom_syl = 'y';
        bool_apf = 'n'; % Does not work with APF
        bool_sntln = 'n'; % Does not work with SNTLN
        fprintf('SNTLN and APF only work when including Matrix Q in sylvester matrix.\n')
        fprintf('Denominator must be included when Q is included \n')
        
end

PrintGlobalVariables();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Add subdirectories

addpath 'Examples'
addpath 'Root Finding Methods'
addpath 'BernsteinMethods'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problemType
    case 'fromRoots'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Get exact polynomial roots from the example file.
        [f_roots_exact] = Examples_Roots(ex_num);
        exact_roots = sortrows(f_roots_exact,1);
        
        % Print the exact roots and multiplicities to terminal
        fprintf('\nExact Roots of Input Polynomial \n');
        fprintf('\t \t \t \t \t Root \t \t \t \t\t  Multiplicity \n')
        fprintf('%30.15f %30.15f \t \t\n',[exact_roots(:,1),exact_roots(:,2)]');
        fprintf('\n');
        
        f_exact_bi          = B_poly(f_roots_exact);
        
        % Get degree of f
        m = length(f_exact_bi) - 1;
        
        % Display the degree of the input polynomial
        disp('Degree of Input Polynomial F ');
        disp(int2str(m));
        
        % Get the Binomial coefficients corresponding to the coefficients of
        % polynomial f.
        Bi_m = zeros(m+1,1);
        for i=1:1:m+1
            Bi_m(i) = nchoosek(m,i-1);
        end
        
        % Get coefficients in Bernstein Basis
        f_exact = f_exact_bi./Bi_m;
        
    case 'fromCoefficients'
        switch ex_num
            case '1'
                f_exact = ...
                    [
                    -0.9865
                    2.2398
                    2.8950
                    1.9092
                    -0.1477
                    ];
            otherwise
                error('Not a valid example number for the *from coefficients* examples.')
        end
end

%%
% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = VariableNoise(f_exact,emin,emax,seed);

%%

% This section calculates the roots by several methods


% Calculate roots by mymethod.
clc_roots_mymthd = o_roots_mymethod(fx);

% Calculate roots by matlab 'roots' function.
clc_roots_mtlb = o_roots_matlab(fx);

% Calculate roots by 'multroot' function.
clc_roots_mltrt = o_roots_multroot(fx);

% Calculate roots by 'Interval Bisection' function
%clc_roots_intvlBsctn = o_roots_bisection(fx);

% Calculate roots by 'Subdivisiton' Method
%clc_roots_subdivision = o_roots_subdivision(ex_num,emin,emax,seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%

% Get vector of exact roots, and extract multiplicities.
% eg: if root one has multiplicity 5, then the vector would be given by
% X1 = [r1 r1 r1 r1 r1 ...].
switch problemType
    case 'fromRoots'
        % Let sum_mult_exct be the sum of all multiplicities of the exact roots
        sum_mult_exct = sum(exact_roots(:,2));
        
        % Initialise a vector to store the roots where roots with high multiplicity
        % are repeated.
        nondistinctRoots_exct = zeros(sum_mult_exct,1);
        
        % Initialise a count
        count = 1;
        
        % for each exact root r_{i}
        for i = 1:1:size(exact_roots,1)
            
            % Get multiplicty of root i
            m = exact_roots(i,2);
            
            % for j = 1,...,m
            for j = 1:1:m
                
                % Add the root to a vector of nondistinct roots
                nondistinctRoots_exct(count,1) = exact_roots(i,1);
                
                % Increment the counter
                count = count + 1;
            end
        end
    case 'fromCoefficients'
    otherwise
        error('Problem type is either fromRoots or fromCoefficients')
end
% Get vector of roots calculated by my method, and extract multiplicities.
% eg: if r_{1} has multiplicity 5, then the vector would be given by
% X1 = [r1 r1 r1 r1 r1 ...].


nondistinctRoots_mymthd = GetRepeatedRoots(clc_roots_mymthd);
nondistinctRoots_mtlb   = GetRepeatedRoots(clc_roots_mtlb);
nondistinctRoots_mltrt  = GetRepeatedRoots(clc_roots_mltrt);



%
% Get the roots by interval bisection
%
try
    % Let sum_root_mult_bsctn be the sum of all of the multiplicities of all of
    % the roots obtained by my bisection method
    sum_root_mult_bsctn = sum(clc_roots_intvlBsctn(:,2));
    
    % Initialise a vector to store the nondistinct roots
    nondistinctRoots_bsctn = zeros(sum_root_mult_bsctn,1);
    
    % Initialise a counter
    count = 1;
    
    for i = 1:1:size(clc_roots_intvlBsctn,1)
        
        % Get multiplicty of root i
        m = clc_roots_intvlBsctn(i,2);
        
        % for each of the m roots at r_{i}
        for j = 1:1:m
            
            % Ad the root to a vector of nondistinct roots
            nondistinctRoots_bsctn(count,1) = clc_roots_intvlBsctn(i,1);
            
            % Increment the counter
            count = count + 1;
        end
    end
catch
end


%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
PLOT_GRAPHS = 'y';
switch PLOT_GRAPHS
    case 'y'
        figure('name','Plot Calculated Roots')
        switch problemType
            case 'fromRoots'
                scatter(real(nondistinctRoots_exct),imag(nondistinctRoots_exct),'black','*','LineWidth',20,'DisplayName','Exact Roots');
        end
        hold on;
    
        scatter((real(nondistinctRoots_mymthd)),imag(nondistinctRoots_mymthd),'yellow','*','DisplayName','My Method');
        scatter((real(nondistinctRoots_mtlb)),imag(nondistinctRoots_mtlb),'red','DisplayName','Matlab Roots');
        scatter((real(nondistinctRoots_mltrt)),imag(nondistinctRoots_mltrt),'green','s','filled','DisplayName','MultRoots');
        xlabel('Real');
        ylabel('Imaginary');
        legend(gca,'show')
        str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
        title(str);
        hold off
    case 'n'
        % Dont plot graph
    otherwise
        error('error: plot_graphs is either y or n')
end

% #########################################################################

% Error Measure One.
% Take the roots, and form a polynomial from the roots calculated by each
% corresponding method. then take the difference in the coefficients of the
% exact polynomial - calculated polynomial / exact polynomial
% \hat{a}_{i} - a_{i} ./ \hat{a}_{i}

% #########################################################################

% Having obtained my root estimates. build the polynomial from the
% calculated roots of the methods above.
calculated_poly_mymthd = B_poly(clc_roots_mymthd);
calculated_poly_mymthd = calculated_poly_mymthd./calculated_poly_mymthd(1) ;


f_exact_bi = B_poly(f_roots_exact);
f_exact_bi  = f_exact_bi ./ f_exact_bi(1);

switch 'problemType'
    case 'fromRoots'
        % Get error measure
        err_mymthd  = (calculated_poly_mymthd - f_exact_bi) ./ f_exact_bi;
        fprintf('\nObtaining polynomial coefficients from calculated roots...\n');
        fprintf('Normalised error in coefficients\n\n')
        
        fprintf('My Method: %g \n\n',norm(err_mymthd))
        
        % Having obtained MultRoot Estimates, build the polynomial from the calculated roots.
        
        calculated_poly_multroot = B_poly(clc_roots_mltrt);
        calculated_poly_multroot = calculated_poly_multroot./calculated_poly_multroot(1);
        
        err_mltrt = ((calculated_poly_multroot - f_exact_bi)) ./ f_exact_bi;
        fprintf('Multroot Method: %g \n\n',norm(err_mltrt));
        
        % Having obtained the Matlab ROOTS, build the polynomial from the calculated roots
        
        calculated_poly_matlabroot = B_poly(clc_roots_mtlb);
        calculated_poly_matlabroot = calculated_poly_matlabroot ./ calculated_poly_matlabroot(1);
        
        err_mtlbrt = ((calculated_poly_matlabroot) - f_exact_bi) ./ f_exact_bi;
        fprintf('MATLAB Roots Method: %g \n\n', norm(err_mtlbrt));
    case 'fromCoefficients'
end




end


function [nondistinctRoots_mymthd] = GetRepeatedRoots(clc_roots_mymthd) 

% Let sum_rt_mult_mymthd be the sum of all of the multiplicities of all of the
% roots obtained by my method
sum_rt_mult_mymthd = sum(clc_roots_mymthd(:,2));

% Initialise a vector to store the nondistinct roots
nondistinctRoots_mymthd = zeros(sum_rt_mult_mymthd,1);

% Initialise a count
count = 1;

% for each unique root i
for i = 1:1:size(clc_roots_mymthd,1)
    
    % Get multiplicty of root i
    m = clc_roots_mymthd(i,2);
    
    % for each of the m roots at r_{i}
    for j = 1:1:m
        
        % Add the root to a vector of nondistinct roots
        nondistinctRoots_mymthd(count,1) = clc_roots_mymthd(i,1);
        
        % Increment the counter
        count = count + 1;
    end
end
end
