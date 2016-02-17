function [] = o_gcd(ex,emin,emax,BOOL_SNTLN,BOOL_APF,BOOL_PREPROC,seed)
% Obtain the Greatest Common Divisor (GCD) of two polynomials defined in
% the example file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs.
%
% ex:   Example Number
%
% emin: Signal to noise ratio (minimum)
%
% emax: Signal to noise ratio (maximum)
%
% BOOL_SNTLN
%
% BOOL_APF
%
% BOOL_DENOM
%
% BOOL_PREPROC
%
% seed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath 'BernsteinMethods'
addpath 'Bezoutian'

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
global bool_sylvesterBuildMethod
global Bool_APFBuildMethod
global geometricMeanMethod
global SEED

bool_bezout = 0;
bool_sntln = BOOL_SNTLN;
bool_apf = BOOL_APF;

bool_denom_syl = 'y';
bool_denom_apf = 'y';

bool_preproc = BOOL_PREPROC;
bool_q = 'y';
bool_log = 'n';

PLOT_GRAPHS = 'n';

MAX_ITERATIONS_SNTLN = 50;
MAX_ITERATIONS_APF = 50;

MAX_ERROR_APF = 1e-10;
MAX_ERROR_SNTLN = 1e-10;

SEED = seed;

bool_sylvesterBuildMethod = 'standard';
Bool_APFBuildMethod = 'standard';
geometricMeanMethod = 'matlab';

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

% output_format (bool)
% the format of output from file o1.m
%   1 - output u v and d in terms of w (coefficients include theta)
%   0 - output u v and d in terms of x
global output_format
output_format = 'dx';

%%
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
%%


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

%% Print the parameters.
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
[f_roots, g_roots,d_roots,t_exact,u_roots,v_roots] = Examples_GCD(ex);

PrintRoots('f',f_roots)
PrintRoots('g',g_roots)
PrintRoots('d',d_roots)

% given the roots of f and g, plot them on a line
switch PLOT_GRAPHS
    case 'n'
        % Dont plot graphs
    case 'y'
        figure('name','Exact roots of f(x), g(x) and d(x)')
        hold on
        title('Roots of f and g on the real interval')
        scatter(f_roots(:,1),ones(size(f_roots(:,1))),'s','DisplayName','Roots of f(x)')
        try
            scatter(g_roots(:,1),ones(size(g_roots(:,1))),'x','DisplayName','Roots of g(x)')
        catch
            fprintf('could not plot exact roots of g\n')
        end
        try
            scatter(d_roots(:,1),ones(size(d_roots(:,1))),'o','DisplayName','Roots of d(x)')
        catch
            fprintf('Could not plot exact roots of d.\n')
        end
        xlabel('Real')
        legend(gca,'show')
        hold off
    otherwise
        error('error PLOT_GRAPH is either y or n')
end


% Display the exact, expected result for the degree of the GCD
fprintf('Degree of GCD of exact input polynomials: %i \n',t_exact)
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
Bi_m = GetBinomials(m);
Bi_n = GetBinomials(n);
Bi_t = GetBinomials(t);
Bi_mt = GetBinomials(m-t);
Bi_nt = GetBinomials(n-t);

% Get exact coefficients of a_{i},b_{i},u_{i},v_{i} and d_{i} of
% polynomials f, g, u, v and d in standard bernstein form.

f_exact = f_exact_bi./Bi_m;
g_exact = g_exact_bi./Bi_n;
d_exact = d_exact_bi./Bi_t;
u_exact = u_exact_bi./Bi_mt;
v_exact = v_exact_bi./Bi_nt;

PrintPoly(f_exact,'f')
PrintPoly(g_exact,'g')

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx = VariableNoise(f_exact,emin,emax,seed);
gx = VariableNoise(g_exact,emin,emax,seed);

% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[~,~,d_calc,u_calc,v_calc] = o1(fx,gx);

% % Normalising the exact values of the gcd, and quotient polynomials.

% Normalise gcd
d_calc = normalise(d_calc);
d_exact = normalise(d_exact);

% Normalise quotient polynomial u
u_calc  = normalise(u_calc);
u_exact = normalise(u_exact);

% Normalise quotient polynomial v
v_calc = normalise(v_calc);
v_exact = normalise(v_exact);


% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.
PrintCoefficients('u',u_exact, u_calc)
PrintCoefficients('v',v_exact, v_calc)
PrintCoefficients('d',d_exact, d_calc)


getError('u',u_calc,u_exact)
getError('v',v_calc,v_exact)
getError('d',d_calc,d_exact)


% Print Errors
fprintf('\nNormwise relative Error in Coefficients \nGiven by: Calculated - exact / exact \n\n');
fprintf('\nNormwise Error in Coefficients  \nGiven by: Calculated - exact \n\n');

end


function [] = PrintRoots(f,f_roots)

% print out the exact roots of f,g and d
fprintf('\nRoots of %s \n',f);
fprintf('\t Root \t \t \t \t\t \t \t   Multiplicity \n')
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',[f_roots(:,1),f_roots(:,2)]');
fprintf('\n');

end

function [] = PrintCoefficients(u,u_exact,u_calc)

fprintf('\nCoefficients of %s \n\n',u);
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')
mat = [real(u_exact(:,1))';  real(u_calc(:,1))' ];
fprintf('%30.15f \t \t \t %30.15f   \t \t \n', mat);
fprintf('\n');

end

function [] = getError(u,u_calc,u_exact)

% Calculate relative errors in outputs
rel_error_uk = norm(abs(u_calc - u_exact) ./ u_exact);

fprintf('\tCalculated relative error %s : %8.2e \n ',u,rel_error_uk);

error_uk = norm(abs(u_calc - u_exact) );

fprintf('\tCalculated error %s : %8.2e \n', u,error_uk);

end

function f = normalise(f)
 
f = f./f(1);

end
