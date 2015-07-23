function [] = o_gcd(ex,emin,emax,BOOL_SNTLN,BOOL_APF,BOOL_DENOM,BOOL_PREPROC,seed)
% Obtain the gcd of two polynomials defined in the example file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ex - (Int) Example Number
% emin - Signal to noise ratio (minimum)
% emax - Signal to noise ratio (maximum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath 'BernsteinMethods'
addpath 'Bezoutian'



global fignum
fignum = 50;


% Initialise global variables to be used across the code

% bool_bezout (boolean)
% 1 - Include bernstein-bezoutian matrix method for comparison purposes.
% 0 - Dont Include bernstein bezoutian matrix
global bool_bezout
bool_bezout = 1;

% BOOL_SNTLN - (Boolean)
%    1 :- Include Structured Perturbations in Sylvester Matrix S_{k} before
%    calculating quotients u_{k}, v_{k} and GCD d_{k}.
%    0 :- Don't Include Structured Perturbations.
global bool_sntln
bool_sntln = BOOL_SNTLN;

% BOOL_APF - (Boolean)
%    1 :- Apply Structured Perturbations to the Approximate Polynomial
%    Factorisation, before obtaining GCD d_{k}.
%    0 :- Don't Include Structured Perturbations.
global bool_apf
bool_apf = BOOL_APF;

% BOOL_DENOM - (Boolean) Given the rearrangement of the Sylvester matrix in
% the Bernstein basis, each partition of each subresultant has a common
% divisor to its elements.
%    1 :- Include Common Denominator.
%    0 :- Exclude Common Denominator.
global bool_denom_syl
bool_denom_syl = BOOL_DENOM;

% BOOL_PREPROC - (Boolean) It has been shown that the inclusion of
% preprocessors Geometric mean, scaling by alpha, change of independent
% variable, yield improved results.
%   1 :- Include Preprocessors.
%   0 :- Exclude Preprocessors.
global bool_preproc
bool_preproc = BOOL_PREPROC;

% BOOL_Q - (Boolean) Consists of binomial coefficients of coprime
% polynomials.
%   1 :- Include the binomial coefficients from the null space in the
%    Sylvester Matrix.
%   0 :- Exclude the binomial coefficients from the null space in the
%    Sylvester Matrix.
global bool_q
bool_q = 1;

% BOOL_LOG - (Boolean)
%   1 :- Perform combinatorial calculations by log method
%   0 :- Perform combinatorial calculations by standard method.
global bool_log
bool_log = 1;

% bool_plotgraphs (boolean)
% 1 - Plot graphs for computations of calculating gcd
% 0 - Avoid plotting graphs (Speed up)
global bool_plotgraphs
bool_plotgraphs = 1;

% set maximum error for sntln problem
global max_error
max_error = 1e-15;

% set maximum number of iterations in the SNTLN and APF problem
global max_iterations
max_iterations = 50;

% bool_sylvesterBuildMethod (boolean)
% 1 :   Build based on individual elements of the Sylvester matrix, each
%       (i,j) element is calculated independently.
% 0 :   Use Naive method calculate D, calculate S, calculate Q, then
%       calculate DTQ.
global bool_sylvesterBuildMethod
bool_sylvesterBuildMethod = 1;

% geometric mean method (bool)
% 0 - Matlab Method
% 1 - My Geometric Mean
% 2 - Set Geometric Mean = 1
global geometricMeanMethod
geometricMeanMethod = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nominal value used when:
% Only one sylvester subresultant exists, ie k = 1 and min(m,n) = 1. where
% m and n are the degrees of input polynomials f and g.
% if max_r./min_r > nominal_value (then minimum value is significantly
% small, to assume that the sylvester matrix is rank deficient)
% then degree is one. otherwise degree is zero
global nominal_value
nominal_value = 100;

% let x be the maximum change in ratio_maxmin_rowsum vector if abs(x) <
% nominal_value_2, if the change is minimal, then all subresultants should
% be classed as rank deficient.

global min_delta_mag_rowsum
min_delta_mag_rowsum = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output_format (bool) 
% the format of output from file o1.m
%   1 - output u v and d in terms of w (coefficients include theta)
%   0 - output u v and d in terms of x 
global output_format
output_format = 0;

% Bool_APFBuildMethod
% 1 :   Build based on individual elements of the Sylvester matrix, each
%       (i,j) element is calculated independently.
% 0 :   Use Naive method calculate D, calculate S, calculate Q, then
%       calculate DTQ.
global Bool_APFBuildMethod
Bool_APFBuildMethod = 1;

% bool_SNTLN_ROOTS
%   1 - Use Roots based SNTLN method, which has the added constraints that
%   g is the derivative of f.
%   0 - Use standard SNTLN where f and g are unconstrained
global bool_SNTLN_Roots 
bool_SNTLN_Roots = 0;

% bool_APF_Roots
%   1 - Use roots based APF method, which has added constraings.
%   0 - use standard apf method where f and g are unconstrained.
global bool_APF_Roots
bool_APF_Roots = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                Consistency of input parameters.

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
% applicable, and the common denominators can not be removed.
% Simplest method, no structure added.
% Override users input options if incompatable.
if (bool_q == 0)
    bool_denom_syl = 1;
    bool_apf = 0; % Does not work with code block APF (Addition of structured perturbation code doesnt exist for exclusion of Q from coefficient matrix).
    bool_sntln = 0; % Does not work with code block SNTLN (Addition of structured perturbations code doesnt exist for exclusion of Q from coefficient matrix).
    fprintf('\nSNTLN and APF only work when including Matrix Q in sylvester matrix.\n')
    fprintf('Denominator must be included when excluding matrix Q \n');
end

% Print the parameters.
fprintf('--------------------------------------------------------------------------- \n')
fprintf('PARAMETERS:\n')
fprintf('\n')
fprintf('\tExample Number: %i \n',ex);
fprintf('\tmin noise : %i \n\tmax noise : %i',emin,emax)
fprintf('\n\tSNTLN : %i \n\tAPF : %i \n\tDENOM : %i \n\tPREPROC : %i \n\tLOG: %i\n\tQ : %i\n',bool_sntln,bool_apf,bool_denom_syl,bool_preproc,bool_log,bool_q);
fprintf('--------------------------------------------------------------------------- \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails

% Add neccesary paths.
addpath 'Measures'
addpath 'Examples'


% Get roots from example file
[f_roots, g_roots,d_roots,rankloss,u_roots,v_roots] = PreviousExamples(ex,seed);

% print out the exact roots of f,g and d
fprintf('\nRoots of f \n');
fprintf('\t Root \t \t \t \t\t \t \t   Multiplicity \n')
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',[f_roots(:,1),f_roots(:,2)]');
fprintf('\n');

fprintf('\nRoots of g \n');
fprintf('\t Root \t \t \t \t\t \t \t   Multiplicity \n')
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',[g_roots(:,1),g_roots(:,2)]');
fprintf('\n');

fprintf('\nRoots of d \n');
fprintf('\t Root \t \t \t \t\t \t \t   Multiplicity \n')
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',[d_roots(:,1),d_roots(:,2)]');
fprintf('\n');

vec_dist = zeros(length(f_roots),length(g_roots));

% for each root in f, get the minimum distance to roots in g
for i = 1:1:length(f_roots)
    for j = 1:1:length(g_roots)
        vec_dist(i,j) = f_roots(i,1) - g_roots(j,1) ;
    end
end

% given the roots of f and g, plot them on a line
figure(50)
hold on
title('Roots of f and g on the real interval')
scatter(f_roots(:,1),ones(size(f_roots(:,1))),'s','DisplayName','Roots of f')
scatter(g_roots(:,1),ones(size(g_roots(:,1))),'x','DisplayName','Roots of g')
scatter(d_roots(:,1),ones(size(d_roots(:,1))),'o','DisplayName','Roots of GCD')
xlabel('Real')

legend(gca,'show')
hold off


% Display the exact, expected result for the degree of the GCD
fprintf('Degree of GCD of exact input polynomials: %i \n',rankloss)
fprintf('--------------------------------------------------------------------------- \n')



% Using roots stored as f and g and obtain polys in scaled bernstein basis
% B_poly returns coefficients $a_{i}$\binom{m}{i} in a scaled bernstein basis.
% We deal with bernstein basis so wish to remove the (mchoosei) such that
% we have $a_{i}$ only which is the coefficient in the Bersntein Basis..
%   fx_exact_bi = \hat{a}_{i} binom{m}{i}
%   gx_exact_bi = \hat{b}_{i} binom{n}{i}

f_exact_bi = B_poly(f_roots);
g_exact_bi = B_poly(g_roots);
d_exact_bi = B_poly(d_roots);
u_exact_bi = B_poly(u_roots);
v_exact_bi = B_poly(v_roots);

% Get degree of polynomials f.
m = length(f_exact_bi) -1;

% Get degree of polynomials g.
n = length(g_exact_bi) -1;

% Get degree of exact GCD
t = length(d_exact_bi) - 1;

% Get sets of binomial coefficients corresponding to each vector

% Bi_m - Binomials corresponding to polynomial f.
Bi_m = zeros(m+1,1);

% Bi_n - Binomials corresponding to polynomial g.
Bi_n = zeros(n+1,1);

% Bi_t - Binomials corresponding to polynomial d.
Bi_t = zeros(t+1,1);

% Bi_nt - Binomails corresponding to polynomial v.
Bi_nt = zeros(n-t+1,1);

% Bi_mt - Binomials corresponding to polynomial u.
Bi_mt = zeros(m-t+1,1);

% supress warnings regarding the nchoosek

% Get the binomial coefficients corresponding to the coefficients of
% f,g,d,v,u.
for i=1:1:m+1
    Bi_m(i) = nchoosek(m,i-1);
end
for i=1:1:n+1
    Bi_n(i) = nchoosek(n,i-1);
end
for i=1:1:t+1
    Bi_t(i) = nchoosek(t,i-1);
end
for i=1:1:m-t+1
    Bi_mt(i) = nchoosek(m-t,i-1);
end
for i=1:1:n-t+1
    Bi_nt(i) = nchoosek(n-t,i-1);
end

% Get exact coefficients of a_{i},b_{i},u_{i},v_{i} and d_{i} of
% polynomials f, g, u, v and d in standard bernstein form.

f_exact = f_exact_bi./Bi_m;
g_exact = g_exact_bi./Bi_n;
d_exact = d_exact_bi./Bi_t;
u_exact = u_exact_bi./Bi_mt;
v_exact = v_exact_bi./Bi_nt;



% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
% fx = $\hat{a}_{i}  + delta\hat{a}_{i}  == a_{i}$
% gx = $\hat{b}_{i}  + delta\hat{b}_{i}  == b_{i}$

fx = VariableNoise(f_exact,emin,emax,seed);
gx = VariableNoise(g_exact,emin,emax,seed);

% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[f_calc,g_calc,d_calc,u_calc,v_calc] = o1(fx,gx);


% output_format (bool)
% The format of output 
%   1 - output u v and d in terms of w (coefficients include theta)
%   0 - output u v and d in terms of x
switch output_format
    case 0 % output in terms of fx,gx,dx 
        f_calc = f_calc;
        g_calc = g_calc;
        u_calc = u_calc;
        v_calc = v_calc;
        d_calc = d_calc;
        
    case 1 % output in terms of fw,gw,dw
        f_calc = f_calc ./ (theta.^vecm);
        g_calc = g_calc ./ (theta.^vecn);
        u_calc = u_calc ./ (theta.^vecmk);
        v_calc = v_calc ./ (theta.^vecnk);
        d_calc = d_calc ./ (theta.^vecnk);
        
end

switch bool_bezout
    % Obtain the gcd by bezoutian method
    case 1
        fprintf('Using Bezout Method')
        o_Bezout(fx,gx)
    case 0
        
end

% % Normalising the exact values of the gcd, and quotient polynomials.
% Normalise gcd
d_calc = d_calc./d_calc(1);
d_exact = d_exact./d_exact(1);

% Normalise quotient polynomial u
u_calc = u_calc./u_calc(1)';
u_exact = u_exact./u_exact(1);

% Normalise quotient polynomial v
v_calc = v_calc./v_calc(1)';
v_exact = v_exact./v_exact(1);

% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.

fprintf('\nCoefficients of u \n\n');
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')
mat = [real(u_exact(:,1))';  real(u_calc(:,1))' ];
fprintf('%30.15f \t \t \t %30.15f   \t \t \n', mat);
fprintf('\n');


fprintf('\nCoefficients of v \n\n');
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')
mat = [real(v_exact(:,1))';  real(v_calc(:,1))' ];
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',mat);
fprintf('\n');


fprintf('\nCoefficients of d \n\n');
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')
mat = [real(d_exact(:,1))';  real(d_calc(:,1))' ];
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',mat);
fprintf('\n');



% Calculate relative errors in outputs
rel_error_uk = norm(abs(u_calc - u_exact) ./ u_exact);
rel_error_vk = norm(abs(v_calc - v_exact) ./ v_exact);
rel_error_dk = norm(abs(d_calc - d_exact) ./ d_exact);

% Print Errors
fprintf('\nNormwise relative Error in Coefficients \nGiven by: Calculated - exact / exact \n\n');
fprintf('\tCalculated relative error u : %8.2e \n', rel_error_uk);
fprintf('\tCalculated relative error v : %8.2e \n', rel_error_vk);
fprintf('\tCalculated relative error d : %8.2e \n\n', rel_error_dk);

% Calculate absolute errors
error_uk = norm(abs(u_calc - u_exact) );
error_vk = norm(abs(v_calc - v_exact) );
error_dk = norm(abs(d_calc - d_exact) );

% Print Errors
fprintf('\nNormwise Error in Coefficients  \nGiven by: Calculated - exact \n\n');
fprintf('\tCalculated error u : %8.2e \n', error_uk);
fprintf('\tCalculated error v : %8.2e \n', error_vk);
fprintf('\tCalculated error d : %8.2e \n\n', error_dk);




end





