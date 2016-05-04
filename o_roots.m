function [] = o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
% o_roots(ex_num,emin,emax,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)
%
% Given an example number, and a set of input parameters, calculate the
% roots r_{i} of the polynomial f(x) and the corresponding multiplicities 
% m_{i}. 
%
% % Inputs
%
% ex_num : (string) Example Number
%
% emin : (float) Noise/Signal maximum threshold (minimum)
%
% emax : (float) Noise/Signal maximum threshold (maximum)
%
% mean_method : (string) method used to compute the mean of entries in C_{n-k}(f)
%               and C_{m-k}(g)
%               'None' - No mean
%               'Geometric Mean Matlab Method'
%               'Geometric Mean My Method'
%              
%
% bool_alpha_theta : (string) {'y,'n''}
%
% low_rank_approx_method : (string) {'None','Standard STLN', 'Standard SNTLN'}
%
% apf_method ('string') {'None', 'Standard APF', 'Root Specific APF'}
%
% % Example
%
% >> o_roots('1',1e-12,1e-10,'Geometric Mean Matlab Method','y','Standard SNTLN','Standard APF')



SetGlobalVariables(ex_num,mean_method,bool_alpha_theta,low_rank_approx_method,apf_method)

global problemType 
problemType = 'fromRoots'; % fromRoots/fromCoefficients


% Validate Inputs.


% Check that max and min signal to noise ratio are the correct way around.
% If not, rearrange min and max.
if emin > emax
    fprintf('minimum noise greater than maximum noise \n swapping values...\n')
    emin_wrong = emin;
    emax_wrong = emax;
    emin = emax_wrong;
    emax = emin_wrong;
end

PrintGlobalVariables();

% Add subdirectories
addpath 'Root Finding Methods'

% Get the polynomial f(x) as a column vector of coefficients.
[fx_exact] = Examples_Roots(ex_num);

% Add Noise to coefficients of exact polynomial f_exact, to obtain noisy
% polynomial fx.
fx = VariableNoise(fx_exact,emin,emax);

% Calculate roots by mymethod.
clc_roots_mymethod = o_roots_mymethod(fx);

% Calculate roots by matlab 'roots' function.
clc_roots_matlab = o_roots_matlab(fx);

% Calculate roots by 'multroot' function.
clc_roots_multroot = o_roots_multroot(fx);

% Calculate roots by 'Interval Bisection' function
% clc_roots_intervalBisection = o_roots_bisection(fx);

% Calculate roots by 'Subdivisiton' Method
%clc_roots_subdivision = o_roots_subdivision(fx);

% Calculate roots by bezier clipping
%clc_roots_clipping = o_roots_BezierClipping(fx);

fx_mymethod = GetWithoutBinomials(B_poly(clc_roots_mymethod));
fx_matlab   = GetWithoutBinomials(B_poly(clc_roots_matlab));
fx_multroot = GetWithoutBinomials(B_poly(clc_roots_multroot));

error_mymethod = norm(Normalise(fx_mymethod) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
error_matlab   = norm(Normalise(fx_matlab) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));
error_multroot = norm(Normalise(fx_multroot) - Normalise(fx_exact)) ./ norm(Normalise(fx_exact));


% %


% Get vector of roots calculated by my method, and extract multiplicities.
% eg: if r_{1} has multiplicity 5, then the vector would be given by
% X1 = [r1 r1 r1 r1 r1 ...].


nondistinctRoots_mymethod = GetRepeatedRoots(clc_roots_mymethod);
nondistinctRoots_matlab   = GetRepeatedRoots(clc_roots_matlab);
nondistinctRoots_multroot  = GetRepeatedRoots(clc_roots_multroot);
nondistinctRoots_bisection  = GetRepeatedRoots(clc_roots_intervalBisection);


%
% Plot the graph real (x) imaginary (y) components of the nondistinct roots
% obtained by the root calculating methods.
global PLOT_GRAPHS
switch PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Plot Calculated Roots',mfilename);
        figure('name',figure_name)
        hold on;
    
        scatter((real(nondistinctRoots_mymethod)),imag(nondistinctRoots_mymethod),'yellow','*','DisplayName','My Method');
        scatter((real(nondistinctRoots_matlab)),imag(nondistinctRoots_matlab),'red','DisplayName','Matlab Roots');
        scatter((real(nondistinctRoots_multroot)),imag(nondistinctRoots_multroot),'green','s','filled','DisplayName','MultRoots');
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
calculated_poly_mymthd = B_poly(clc_roots_mymethod);
calculated_poly_mymthd = calculated_poly_mymthd./calculated_poly_mymthd(1) ;


fx_bi = B_poly(f_roots_exact);
fx_bi = Normalise(fx_bi);

switch 'problemType'
    case 'fromRoots'
        % Get error measure
        err_mymthd  = (calculated_poly_mymthd - fx_bi) ./ fx_bi;
        fprintf('\nObtaining polynomial coefficients from calculated roots...\n');
        fprintf('Normalised error in coefficients\n\n')
        
        fprintf('My Method: %g \n\n',norm(err_mymthd))
        
        % Having obtained MultRoot Estimates, build the polynomial from the calculated roots.
        
        calculated_poly_multroot = B_poly(clc_roots_multroot);
        calculated_poly_multroot = calculated_poly_multroot./calculated_poly_multroot(1);
        
        err_mltrt = ((calculated_poly_multroot - fx_bi)) ./ fx_bi;
        fprintf('Multroot Method: %g \n\n',norm(err_mltrt));
        
        % Having obtained the Matlab ROOTS, build the polynomial from the calculated roots
        
        calculated_poly_matlabroot = B_poly(clc_roots_matlab);
        calculated_poly_matlabroot = calculated_poly_matlabroot ./ calculated_poly_matlabroot(1);
        
        err_mtlbrt = ((calculated_poly_matlabroot) - fx_bi) ./ fx_bi;
        fprintf('MATLAB Roots Method: %g \n\n', norm(err_mtlbrt));
    case 'fromCoefficients'
end

% % 
PrintToFile(m);


end


function [nondistinctRoots_mymthd] = GetRepeatedRoots(mat_Root_Mult) 
% Given the matrix whose columns are a [root,multiplicity] pair, get a
% vector which contains each root r_{i} m_{i} times, where m_{i} is the
% multiplicity of r_{i}.


% Let sum_rt_mult_mymthd be the sum of all of the multiplicities of all of the
% roots obtained by my method
sum_rt_mult_mymthd = sum(mat_Root_Mult(:,2));

% Initialise a vector to store the nondistinct roots
nondistinctRoots_mymthd = zeros(sum_rt_mult_mymthd,1);

% Initialise a count
count = 1;

% for each unique root i
for i = 1:1:size(mat_Root_Mult,1)
    
    % Get multiplicty of root i
    m = mat_Root_Mult(i,2);
    
    % for each of the m roots at r_{i}
    for j = 1:1:m
        
        % Add the root to a vector of nondistinct roots
        nondistinctRoots_mymthd(count,1) = mat_Root_Mult(i,1);
        
        % Increment the counter
        count = count + 1;
    end
end
end


function []= PrintToFile(m)

global NOISE
global BOOL_PREPROC
global LOW_RANK_APPROXIMATION_METHOD
global APF_METHOD

fullFileName = 'o_roots_results.txt';


if exist(fullFileName, 'file')
    fileID = fopen(fullFileName,'a');
    fprintf(fileID,'%5d \t %s \t %5s \t %s \t %s \n',...
        m,NOISE,BOOL_PREPROC, LOW_RANK_APPROXIMATION_METHOD,APF_METHOD);
    fclose(fileID);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
  uiwait(msgbox(warningMessage));
end




end
