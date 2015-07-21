function [] = o_roots(ex,emin,emax,BOOL_SNTLN,BOOL_APF,BOOL_DENOM,BOOL_PREPROC,seed)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% ex - (Int) Example Number

% emin - Noise/Signal maximum threshold (minimum)

% emax - Noise/Signal maximum threshold (maximum)

% BOOL_SNTLN - Assigned to global variable (see below)

% BOOL_APF - Assigned to global variable (see below)

% BOOL_DENOM - Assigned to global variable (see below)

% Bool_Preproc - Assigned to global variable (see below)

% seed (int) - Integer chosen to randomly generate the roots and
% multiplicities of input polynomials f and g as well as the noise which is
% added to their coefficients.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                    # Global Variables #

% BOOL_SNTLN - (Boolean)
%    1 :- Include Structured Perturbations in Sylvester Matrix S_{k} before
%    calculating quotients u_{k}, v_{k} and GCD d_{k}.
%    0 :- Don't Include Structured Perturbations.

% BOOL_APF - (Boolean)
%    1 :- Apply Structured Perturbations to the Approximate Polynomial
%    Factorisation, before obtaining GCD d_{k}.
%    0 :- Don't Include Structured Perturbations.

% bool_denom_apf - (Boolean)
% Given the rearrangement of the coefficient matrix in the approximate
% polynomial factorisation [C(u);C(v)]d = [f;g]
% The entries in each partition C(u) and C(v) have a common divisor, which
% is removed by manipulation of the binomial coefficients.
% NOTE - This should always be included.
%    1 :- Include Common Denominators.
%    0 :- Exclude Common Denominators.

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.

% BOOL_PREPROC - (Boolean) It has been shown that the inclusion of
% preprocessors Geometric mean, scaling by alpha, change of independent
% variable, yield improved results.
%   1 :- Include Preprocessors.
%   0 :- Exclude Preprocessors.
global bool_sntln
global bool_apf
global bool_denom_apf
global bool_denom_syl
global bool_preproc

bool_sntln = BOOL_SNTLN;
bool_apf = BOOL_APF;
bool_denom_apf = 1;
bool_denom_syl = BOOL_DENOM;
bool_preproc = BOOL_PREPROC;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                            CONSTANTS
global fignum
fignum = 50;

% bool_SNTLN_ROOTS
% 1 :   Use Roots based SNTLN method, which has the added constraints that
%       g is the derivative of f.
% 0 :   Use standard SNTLN where f and g are unconstrained

% bool_APF_Roots
% 1 :   Use roots based APF method, which has added constraings.
% 0 :   Use standard apf method where f and g are unconstrained.

% bool_q - (Boolean) Consists of binomial coefficients of coprime
% polynomials.
% 1 :   Include the binomial coefficients from the null space in the
%       Sylvester Matrix.
% 0 :   Exclude the binomial coefficients from the null space in the
%       Sylvester Matrix.

% BOOL_LOG - (Boolean)
% 1 :   Perform calculations by log method
% 0 :   Perform calculations by standard method.

% BOOL_deconvolve - (Boolean)
% 0 - Standard deconvolution in the root finder
% 1 - Batch deconvolution in the root finder

% max_error - set the maximum error allowed in the two LSE Problems, the
%   first of which
%   is the SNTLN of the Sylvester matrix, the second is the approximate
%   polynomial factorisation.

% max_iterations - Set the maximum number of iterations in the two LSE
%   problems.

% bool_reordercols (bool)
% 0 :   Leave columns of the Sylvester Matrix as standard partitions, where
%       the first n-k+1 columns contain entries corresponding to the
%       coefficients of f(y), and the last m-k+1 columns correspond to
%       coefficients of g(y)
% 1 :   Rearrange columns of the Sylvester subresultant matrices in
%       accordance with Z Zeng - Computing Multiple Roots of inexact
%       polynomials (page 889)

% plotgraphs (bool)
% 0 :   Don't plot graphs, just perform root finding operation.
% 1 :   Plot Graphs associated with calculating the GCD

% Bool_sylvesterBuildMethod
% 1 :   Build based on individual elements of the Sylvester matrix, each
%       (i,j) element is calculated independently.
% 0 :   Use Naive method calculate D, calculate S, calculate Q, then
%       calculate DTQ.

% Bool_APFBuildMethod
% 1 :   Build based on individual elements of the Sylvester matrix, each
%       (i,j) element is calculated independently.
% 0 :   Use Naive method calculate D, calculate S, calculate Q, then
%       calculate DTQ.

global bool_SNTLN_Roots
global bool_APF_Roots
global bool_q
global bool_log
global bool_deconvolve
global max_error
global max_iterations
global bool_plotgraphs
global bool_sylvesterBuildMethod
global Bool_APFBuildMethod

bool_SNTLN_Roots = 1;
bool_APF_Roots = 1;
bool_q = 1;
bool_log = 1;
bool_deconvolve = 1;
max_error = 1e-15;
max_iterations = 50;
bool_plotgraphs = 0;

bool_sylvesterBuildMethod = 1;
Bool_APFBuildMethod = 1;

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
min_delta_mag_rowsum = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output_format (bool)
% The format of output from file o1.m
%   1 - output u v and d in terms of w (coefficients include theta)
%   0 - output u v and d in terms of x
global output_format
output_format = 1;


% geometricMeanMethod
% used when calculating geometric means of the entries of a Sylvester
% matrix, a standard method would consider each entry in the matrix, but a
% new method, described in my internal report offers speed up due to the
% structured nature of the Sylvester matrix entries for the Sylvester
% matrix in the Bernstein basis.
% 0 -   use MatLab Built in method for calculating geometric means
% 1 -   use my method of calculating geometric means#
% 2 -   Set geometric mean equal to one, essentially eliminating the
%       preprocessor.
global geometricMeanMethod
geometricMeanMethod = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if (bool_q == 0)
    
    bool_denom_syl = 1;
    bool_apf = 0; % Does not work with APF
    bool_sntln = 0; % Does not work with SNTLN
    fprintf('SNTLN and APF only work when including Matrix Q in sylvester matrix.\n')
    fprintf('Denominator must be included when Q is included \n')
end

fprintf('PARAMETERS:\n\n')
fprintf('\tExample Number: %i\n',ex)
fprintf('\tmin noise : %i \n\tmax noise : %i\n',emin,emax)
fprintf('INPUT VARIABLES')
fprintf('\n\tSNTLN : %i \n',...
    bool_sntln);
fprintf('\tAPF : %i \n ',bool_apf)
fprintf('\tDENOM : %i \n',bool_denom_syl)
fprintf('\tPREPROC : %i \n',bool_preproc)
fprintf('\tLOG: %i\n',bool_log)
fprintf('\tQ : %i\n',bool_q)
fprintf('\tSNTLN Derivative Constraint: %i \n',bool_SNTLN_Roots)
fprintf('')
fprintf('--------------------------------------------------------------------------- \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Add subdirectories

addpath 'Examples'
addpath 'Root Finding Methods'
addpath 'BernsteinMethods'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get exact polynomial roots from the example file.
[f_roots_exact] = Root_Examples(ex,seed);
exact_roots = sortrows(f_roots_exact,1);

% Print the exact roots and multiplicities to terminal
fprintf('\nExact Roots of Input Polynomial \n');
fprintf('\t \t \t \t \t Root \t \t \t \t\t  Multiplicity \n')
fprintf('%30.15f %30.15f \t \t\n',[exact_roots(:,1),exact_roots(:,2)]');
fprintf('\n');

% Calculate roots by mymethod.
clc_roots_mymthd = o_roots_mymethod(ex,emin,emax,seed);

% Calculate roots by matlab 'roots' function.
clc_roots_mtlb = o_roots_matlab(ex,emin,emax,seed);

% Calculate roots by 'multroot' function.
clc_roots_mltrt = o_roots_multroot(ex,emin,emax,seed);

% Calculate roots by 'Interval Bisection' function
clc_roots_intvlBsctn = o_roots_bisection(ex,emin,emax,seed);

% Calculate roots by 'Subdivisiton' Method
clc_roots_subdivision = o_roots_subdivision(ex,emin,emax,seed);


% Get vector of exact roots, and extract multiplicities.
% eg: if root one has multiplicity 5, then X1 = [r1 r1 r1 r1 r1 ...]

% Let n1 be the total number of non distinct roots
n1 = sum(exact_roots(:,2));

% Initialise a vector to store the nondistinct roots (contains repetition)
nondistinctRoots_exct = zeros(n1,1);

count = 1;
for i = 1:1:size(exact_roots,1)
    % Get multiplicty of root i
    m = exact_roots(i,2);
    for j = 1:1:m
        % Repeat root m times, where m is multplicity
        nondistinctRoots_exct(count,1) = exact_roots(i,1);
        count = count + 1;
    end
end

% Get vector of my calculated roots
% Form a vector of roots

% Let n2 be the total number of nondistinct roots obtained from my method
num_nondistinct_rt_mymthd = sum(clc_roots_mymthd(:,2));

% Initialise a vector to store the nondistinct roots
nondistinctRoots_mymthd = zeros(num_nondistinct_rt_mymthd,1);

% Initialise a count
count = 1;

% for each unique root i
for i = 1:1:size(clc_roots_mymthd,1)
    % Get multiplicty of root i
    m = clc_roots_mymthd(i,2);
    
    % for each of the m roots at i
    for j = 1:1:m
        nondistinctRoots_mymthd(count,1) = clc_roots_mymthd(i,1);
        count = count + 1;
    end
end


% Get vector of roots by matlab method

% Let n3 be the total number of nondistinct roots obtained from matlab
% method
n_nondistinctRoots_mtlb = sum(clc_roots_mtlb(:,2));
nondistinctRoots_mtlb = zeros(n_nondistinctRoots_mtlb,1);
count = 1;
for i = 1:1:length(clc_roots_mtlb)
    % Get multiplicty of root i
    m = clc_roots_mtlb(i,2);
    for j = 1:1:m
        nondistinctRoots_mtlb(count,1) = clc_roots_mtlb(i,1);
        count = count + 1;
    end
end

%
% Get the roots by MULTROOT zheng function
% Form a vector of roots from MULTROOT.
n_nondistinctRoots_mltrt = sum(clc_roots_mltrt(:,2));
nondistinctRoots_mltrt = zeros(n_nondistinctRoots_mltrt,1);
count = 1;
for i = 1:1:size(clc_roots_mltrt,1)
    % Get multiplicty of root i
    m = clc_roots_mltrt(i,2);
    for j = 1:1:m
        nondistinctRoots_mltrt(count,1) = clc_roots_mltrt(i,1);
        count = count + 1;
    end
end

%
% Get the roots by interval bisection
%
try
    n_nondistinctRoots_bsctn = sum(clc_roots_intvlBsctn(:,2));
    nondistinctRoots_bsctn = zeros(n_nondistinctRoots_bsctn,1);
    count = 1;
    
    for i = 1:1:size(clc_roots_intvlBsctn,1)
        % Get multiplicty of root i
        m = clc_roots_intvlBsctn(i,2);
        for j = 1:1:m
            nondistinctRoots_bsctn(count,1) = clc_roots_intvlBsctn(i,1);
            count = count + 1;
        end
    end
catch
end
%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.


figure()
scatter(real(nondistinctRoots_exct),imag(nondistinctRoots_exct),'black*','DisplayName','Exact Roots');
hold on;
scatter(real(nondistinctRoots_mymthd),imag(nondistinctRoots_mymthd),'blueo','DisplayName','My Method');
scatter(real(nondistinctRoots_mtlb),imag(nondistinctRoots_mtlb),'red','DisplayName','Matlab Roots');
scatter(real(nondistinctRoots_mltrt),imag(nondistinctRoots_mltrt),'green*','DisplayName','MultRoots');
xlabel('Real');
ylabel('Imaginary');
str = sprintf('Plot of Calculated Roots of Polynomial f(y). \n componentwise noise = %g',emin);
title(str);
hleg = legend();
hold off


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

[f_roots_exact,~] = Root_Examples(ex,seed);

f_exact_bi = B_poly(f_roots_exact);
f_exact_bi  = f_exact_bi ./ f_exact_bi(1);

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


end

