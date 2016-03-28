function [] = o_roots(ex_num,emin,emax,bool_preproc,low_rank_approx_method,apf_method)
% Given an example number, and a set of input parameters, calculate the
% roots r_{i} of the polynomial f(x) and the corresponding multiplicities 
% m_{i}. 
%
%                           Inputs
%
% ex - (Int) Example Number
%
% emin - Noise/Signal maximum threshold (minimum)
%
% emax - Noise/Signal maximum threshold (maximum)
%
% low_rank_approx_method - Assigned to global variable (see below)
%
% apf_method - {'None', 'Standard APF', 'Root Specific APF'}
%
% bool_preproc - Assigned to global variable (see below)
%
%

SetGlobalVariables(bool_preproc,low_rank_approx_method,apf_method)

global problemType 
problemType = 'fromRoots'; % fromRoots/fromCoefficients


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
        
        f_bi = B_poly(f_roots_exact);
        f_exact = GetWithoutBinomials(f_bi); 
        
        % Get degree of f
        m = length(f_bi) - 1;
        
        % Display the degree of the input polynomial
        disp('Degree of Input Polynomial F ');
        disp(int2str(m));
        
        
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
fx = VariableNoise(f_exact,emin,emax);

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


f_bi = B_poly(f_roots_exact);
f_bi  = f_bi ./ f_bi(1);

switch 'problemType'
    case 'fromRoots'
        % Get error measure
        err_mymthd  = (calculated_poly_mymthd - f_bi) ./ f_bi;
        fprintf('\nObtaining polynomial coefficients from calculated roots...\n');
        fprintf('Normalised error in coefficients\n\n')
        
        fprintf('My Method: %g \n\n',norm(err_mymthd))
        
        % Having obtained MultRoot Estimates, build the polynomial from the calculated roots.
        
        calculated_poly_multroot = B_poly(clc_roots_mltrt);
        calculated_poly_multroot = calculated_poly_multroot./calculated_poly_multroot(1);
        
        err_mltrt = ((calculated_poly_multroot - f_bi)) ./ f_bi;
        fprintf('Multroot Method: %g \n\n',norm(err_mltrt));
        
        % Having obtained the Matlab ROOTS, build the polynomial from the calculated roots
        
        calculated_poly_matlabroot = B_poly(clc_roots_mtlb);
        calculated_poly_matlabroot = calculated_poly_matlabroot ./ calculated_poly_matlabroot(1);
        
        err_mtlbrt = ((calculated_poly_matlabroot) - f_bi) ./ f_bi;
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
