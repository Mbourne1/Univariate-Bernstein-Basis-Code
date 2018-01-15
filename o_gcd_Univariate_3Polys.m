function [] = o_gcd_Univariate_3Polys(ex_num_var, emin, emax, mean_method, ...
    bool_alpha_theta, low_rank_approx_method, apf_method, ...
    sylvester_build_method, nEquations, rank_revealing_metric)
% o_gcd_3Polys(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Obtain the Greatest Common Divisor (GCD) d(x) of two polynomials f(x) and
% g(x) as defined in the example file.
%
%
% % Inputs.
%
% ex : (String) Example Number
%
% emin: (Float) Signal to noise ratio (minimum)
%
% emax: (Float) Signal to noise ratio (maximum)
%
% mean_method : Method for taking mean of entires in S_{k}
%           'Geometric Mean Matlab Method'
%           'Geometric Mean My Method'
%
% bool_alpha_theta : (Bool) if preprocessing is performed
%
% low_rank_approx_method : (String)
%           'Standard STLN'
%           'Standard SNTLN'
%           'Root Specific SNTLN'
%           'None'
%
% apf_method : (String)
%           'Standard APF NonLinear'
%           'Standard APF Linear'
%           'None'
%
% Sylvester_Build_Method : (String)
%           'T'
%           'DT'
%           'DTQ'
%           'TQ'
%           'DTQ Denominator Removed'
%           'DTQ Rearranged'
%
% nEquations : (String)
%
%   '2' : Syvlester matrix defined by 2 equations, and has 2 x 3 structure
%   '3' : Sylvester matrix defined by 3 equations, and has 3 x 3 structure
%
%
%
%
%
%
% % Example
% >> o_gcd_Univariate_3Polys('1a', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'None', 'None', 'DTQ', '2' 'Minimum Singular Values', true)
% >> o_gcd_Univariate_3Polys('1a', 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'Standard APF Nonlinear', 'DTQ', '2', 'Minimum Singular Values', true)
%
% % Custom Example
%
% ex_num = 'Custom:m=10 n=10 t=5 low=0 high=1'
% >> o_gcd_Univariate_3Polys(ex_num, 1e-12, 1e-10, 'Geometric Mean Matlab Method', true, 'Standard STLN', 'None', 'DTQ', 'Minimum Singular Values')


% if not 11 arguments then error
if (nargin ~= 10)
    error('Not enough input arguments')
end

global SETTINGS



% Add subfolders
restoredefaultpath
addpath(genpath(pwd))


% % Ensure that minimum noise level is less than maximum noise level
if emin > emax
    temp = emin;
    emin = emax;
    emax = temp;
end

% % Set global variables.
SetGlobalVariables_GCD_3Polys(...
    ex_num_var,...
    emin,...
    emax,...
    mean_method,...
    bool_alpha_theta,...
    low_rank_approx_method,...
    apf_method,...
    sylvester_build_method, ...
    nEquations,...
    rank_revealing_metric);

% Print the parameters.
LineBreakLarge()
fprintf('INPUT VARIABLES\n');
fprintf('\t EXAMPLE NUMBER : %s \n', ex_num_var);
fprintf('\t EMIN : %s \n' , num2str(SETTINGS.EMIN));
fprintf('\t EMAX : %s \n' , num2str(SETTINGS.EMAX));
fprintf('\t MEAN METHOD : %s \n', SETTINGS.MEAN_METHOD);
fprintf('\t ALPHA_THETA : %s \n', num2str(SETTINGS.BOOL_ALPHA_THETA));
fprintf('\t LOW RANK APPROX METHOD : %s \n', SETTINGS.LOW_RANK_APPROXIMATION_METHOD);
fprintf('\t APF METHOD : %s \n ', SETTINGS.APF_METHOD);
fprintf('\t LOG: %s \n', SETTINGS.BOOL_LOG);
fprintf('\t SYLVESTER BUILD METHOD: %s \n', SETTINGS.SYLVESTER_BUILD_METHOD);
fprintf('\t SYLVESTER n EQUATIONS: %s \n', SETTINGS.SYLVESTER_EQUATIONS);
fprintf('\t SCALING METHOD: %s \n', SETTINGS.SCALING_METHOD);
LineBreakLarge()

% o - gcd - Calculate GCD of two Arbitrary polynomials
% Given two sets of polynomial roots, form polynomials f and g, expressed
% in the Bernstein Basis. Add noise, and calculate the GCD of the two
% polynomails



% Get roots from example file
[fx_exact, gx_exact, hx_exact, dx_exact, ux_exact, vx_exact, wx_exact] = ...
    Examples_GCD_3Polys(ex_num_var);

% Add componentwise noise to coefficients of polynomials in 'Standard Bernstein Basis'.
fx = AddVariableNoiseToPoly(fx_exact, emin, emax);
gx = AddVariableNoiseToPoly(gx_exact, emin, emax);
hx = AddVariableNoiseToPoly(hx_exact, emin, emax);


%PlotCoefficients({fx, gx, hx}, {'fx','gx','hx'});

m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);
t_exact = GetDegree(dx_exact);


LineBreakLarge()
fprintf('\t m : %i \n', m)
fprintf('\t n : %i \n', n)
fprintf('\t o : %i \n', o)
fprintf('\t t : %i \n', t_exact)
LineBreakLarge()


% set upper and lower limit of the degree of the GCD. Since this is a
% general GCD problem, no prior limits are known.
lower_limit_t = 1;
upper_limit_t = min([m n o]);

limits_t = [lower_limit_t, upper_limit_t];
rank_range = [0 0];


% Obtain the coefficients of the GCD d and quotient polynomials u and v.
[~, ~, ~, dx_calc, ux_calc, vx_calc, wx_calc] = ...
    o_gcd_3Polys_mymethod(fx, gx, hx, limits_t, rank_range);

% Check coefficients of calculated polynomials are similar to those of the
% exact polynomials.

LineBreakMedium();

% try
    my_error.dx = GetPolynomialError(dx_exact, dx_calc);
    my_error.ux = GetPolynomialError(ux_exact, ux_calc);
    my_error.vx = GetPolynomialError(vx_exact, vx_calc);
    my_error.wx = GetPolynomialError(wx_exact, wx_calc);
% catch
%     my_error.dx = 1000;
%     my_error.ux = 1000;
%     my_error.vx = 1000;
%     my_error.wx = 1000;
% end

fprintf('Error u(x) : %e \n', my_error.ux);
fprintf('Error v(x) : %e \n', my_error.vx);
fprintf('Error w(x) : %e \n', my_error.wx);
fprintf('Error d(x) : %e \n', my_error.dx);
fprintf('Average Error : %e \n', mean([my_error.ux, my_error.vx, my_error.dx]))


% Print results to results file
%PrintToFile(GetDegree(fx), GetDegree(gx), GetDegree(hx), GetDegree(dx_calc), my_error)
LineBreakMedium();

end




function [] = PrintToFile(m, n, o, t, error)
% Print Results to file
%
% % Inputs
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% o : (Int) Degree of polynomial h(x)
%
% t : (Int) Computed degree of the GCD d(x)
%
% error : array of errors e
%   error.dx : (Float)
%   error.ux
%   error.vx
%   error.wx


global SETTINGS

% Specify file name
fullFileName = sprintf('Results/Results_o_gcd_3Polys.txt');

% If file already exists append a line
if exist(fullFileName, 'file')
    
    fileID = fopen(fullFileName, 'a');
    WriteNewLine()
    fclose(fileID);
    
else % File doesnt exist so create it
    
    fileID = fopen( fullFileName, 'wt' );
    WriteHeader()
    WriteNewLine()
    fclose(fileID);
    
end



    function WriteNewLine()
        
        fprintf(fileID,'%s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s, %s,%s,%s \n',...
            datetime('now'),...
            SETTINGS.EX_NUM,...
            int2str(m),...
            int2str(n),...
            int2str(o),...
            int2str(t),...
            num2str(error.ux),...
            num2str(error.vx),...
            num2str(error.wx),...
            num2str(error.dx),...
            SETTINGS.MEAN_METHOD,...
            num2str(SETTINGS.BOOL_ALPHA_THETA),...
            SETTINGS.EMIN,...
            SETTINGS.EMAX,...
            SETTINGS.LOW_RANK_APPROXIMATION_METHOD,...
            num2str(SETTINGS.LOW_RANK_APPROX_REQ_ITE),...
            SETTINGS.APF_METHOD,...
            num2str(SETTINGS.APF_REQ_ITE),...
            SETTINGS.BOOL_LOG,...
            SETTINGS.SYLVESTER_BUILD_METHOD,...
            SETTINGS.GCD_COEFFICIENT_METHOD...
            );
        % 21 inputs
        
    end

    function WriteHeader()
        fprintf(fileID,'DATE,EX_NUM,m,n,o,t,ERROR_UX,ERROR_VX,ERROW_WX,ERROR_DX,MEAN_METHOD,BOOL_ALPHA_THETA, EMIN, EMAX, LOW_RANK_APPROX_METHOD,LOW_RANK_ITE, APF_METHOD, APF_ITE,BOOL_LOG,SYLVESTER_BUILD_METHOD,GCD_METHOD\n');
    end


end

