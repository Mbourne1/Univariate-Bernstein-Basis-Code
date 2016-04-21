function [] = o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
% o_gcd(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% Inputs.
%
% ex:   Example Number
%
% emin: Signal to noise ratio (minimum)
%
% emax: Signal to noise ratio (maximum)
%
% mean_method : Method for taking mean of entires in S_{k}
%
% bool_alpha_theat : 'y' or 'n' if preprocessing is performed
%
% low_rank_approx_method : 'Standard STLN' 'Standard SNTLN'
%
% apf_method :


addpath 'Bezoutian'

SetGlobalVariables(mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
%                Consistency of input parameters.

% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
if emin > emax
    fprintf('minimum noise greater than maximum noise \n swapping values...\n')
    [emin,emax] = swap(emin,emax);
end



% Print the parameters.
PrintGlobalVariables()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails

% Add neccesary paths.
addpath 'Measures'

% Get roots from example file
[f_roots, g_roots,d_roots,t_exact,u_roots,v_roots] = Examples_GCD(ex_num);

PrintRoots('f',f_roots)
PrintRoots('g',g_roots)
PrintRoots('d',d_roots)

% Plot the roots of f(x), g(x) and d(x)
PlotRoots()

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


% Get exact coefficients of a_{i},b_{i},u_{i},v_{i} and d_{i} of
% polynomials f, g, u, v and d in standard bernstein form.
f_exact = GetWithoutBinomials(f_exact_bi);
g_exact = GetWithoutBinomials(g_exact_bi);
d_exact = GetWithoutBinomials(d_exact_bi);
u_exact = GetWithoutBinomials(u_exact_bi);
v_exact = GetWithoutBinomials(v_exact_bi);

% Print the coefficients of f(x) and g(x)
PrintCoefficients_Bivariate_Bernstein(f_exact,'f')
PrintCoefficients_Bivariate_Bernstein(g_exact,'g')

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx = VariableNoise(f_exact,emin,emax);
gx = VariableNoise(g_exact,emin,emax);

% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[~,~,dx_calc,ux_calc,vx_calc] = o1(fx,gx);

% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.
PrintCoefficients('u',u_exact, ux_calc);
PrintCoefficients('v',v_exact, vx_calc);
PrintCoefficients('d',d_exact, dx_calc);

error_dx = GetError('d',d_exact,dx_calc);

PrintToFile(GetDegree(fx),GetDegree(gx),GetDegree(dx_calc),error_dx)


end


function [] = PrintRoots(f,f_roots)

% print out the exact roots of f,g and d
fprintf('\nRoots of %s \n',f);
fprintf('\t Root \t \t \t \t\t \t \t   Multiplicity \n')
fprintf('%30.15f \t \t \t %30.15f   \t \t \n',[f_roots(:,1),f_roots(:,2)]');
fprintf('\n');

end

function [] = PrintCoefficients(name,u_exact,u_calc)

% Normalise quotient polynomial u
u_calc  = Normalise(u_calc);
u_exact = Normalise(u_exact);

fprintf('\nCoefficients of %s \n\n',name);
fprintf('\t Exact \t \t \t \t\t \t \t   Computed \n')
mat = [real(u_exact(:,1))';  real(u_calc(:,1))' ];
fprintf('%30.15f \t \t \t %30.15f   \t \t \n', mat);
fprintf('\n');
GetError(name,u_calc,u_exact);

end

function [error] = GetError(name,f_calc,f_exact)
% Get distance between f(x) and the calulated f(x)

f_calc  = Normalise(f_calc);
f_exact = Normalise(f_exact);

% Calculate relative errors in outputs
rel_error_f = norm(abs(f_calc - f_exact) ./ f_exact);

fprintf('\tCalculated relative error %s : %8.2e \n ',name,rel_error_f);

error = norm(abs(f_calc - f_exact) );

fprintf('\tCalculated error %s : %8.2e \n', name,error);

end



function []= PrintToFile(m,n,t,error_dx)

global NOISE
global BOOL_PREPROC
global LOW_RANK_APPROXIMATION_METHOD
global APF_METHOD

fullFileName = 'o_gcd_results.txt';


if exist('o_gcd_results.txt', 'file')
    fileID = fopen('o_gcd_results.txt','a');
    fprintf(fileID,'%5d \t %5d \t %5d \t %s \t %s \t %s \t %s \t %s \n',...
        m,n,t,error_dx,BOOL_PREPROC, NOISE, LOW_RANK_APPROXIMATION_METHOD,APF_METHOD);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end




end


