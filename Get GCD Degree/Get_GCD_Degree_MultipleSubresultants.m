
function [t] = Get_GCD_Degree_MultipleSubresultants(vMinimumSingularValues,deg_limits)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.


% Intialise global settings
global SETTINGS

% Get the function which called this function.
[St,~] = dbstack();
calling_function = St(2).name;

lower_lim = deg_limits(1);
upper_lim = deg_limits(2);




% Get the upper bound of the degree of the GCD
min_mn = lower_lim + length(vMinimumSingularValues) - 1;

% Get the maximum change in singular values and the index at which the
% maximum change occured.
[maxChangeSingularValues, indexMaxChange] = Analysis(vMinimumSingularValues);

display([mfilename ' : ' calling_function ' : ' sprintf('Max Change : %2.4f ', maxChangeSingularValues)]);
display([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %2.4f ',SETTINGS.THRESHOLD)]);


if  abs(maxChangeSingularValues) < SETTINGS.THRESHOLD
        
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = log10(mean(vMinimumSingularValues));
    
    % %
    % %
    % %
    switch SETTINGS.PLOT_GRAPHS
        case 'y'
            figure_name = sprintf([mfilename ' : ' calling_function ': Singular Values']);
            figure('name',figure_name)
            plot(log10(vMinimumSingularValues));
            hold on
            mu = avgMinSingularValue;
            hline = refline([0 mu]);
            hline.Color = 'r';  
            vline(lower_lim,'b','');
            vline(lower_lim,'b','');
            hold off
            
        case 'n'
    end
    
    
    if  avgMinSingularValue < SETTINGS.THRESHOLD_RANK
        % If all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        t = min_mn;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Subresultants are Rank Deficient \n')])
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
        
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All subresultants Full Rank \n')])
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
        
    end
else
    % maxChange is signifcant
    t = lower_lim + indexMaxChange - 1;
    
    fprintf([mfilename ' : ' calling_function ' : ' 'Mixed \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
    
    
end

end