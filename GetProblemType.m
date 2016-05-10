
function [t] = GetProblemType(vMinimumSingularValues,lower_lim)
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

% Get the upper bound of the degree of the GCD
min_mn = lower_lim + length(vMinimumSingularValues) - 1;

% Get the maximum change in singular values and the index at which the
% maximum change occured.
[maxChangeSingularValues, indexMaxChange] = Analysis(vMinimumSingularValues);

fprintf([calling_function ' : ' sprintf('Max Change : %2.4f \n', maxChangeSingularValues)]);
fprintf([calling_function ' : ' sprintf('Threshold : %2.4f \n',SETTINGS.THRESHOLD)]);


if  abs(maxChangeSingularValues) < SETTINGS.THRESHOLD
        
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = log10(mean(vMinimumSingularValues));
    
    % %
    % %
    % %
    figure_name = sprintf([calling_function ': Singular Values']);
    figure('name',figure_name)
    plot(log10(vMinimumSingularValues));
    hold on
    mu = avgMinSingularValue;
    hline = refline([0 mu]);
    hline.Color = 'r';  
    hold off
    
    
    if  avgMinSingularValue < SETTINGS.THRESHOLD
        % If all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        fprintf([calling_function ' : ' sprintf('All Subresultants are Rank Deficient \n')])
        t = min_mn;
        
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        fprintf([calling_function ' : ' sprintf('All subresultants Full Rank \n')])
        t = 0;
    end
else
    % maxChange is signifcant
    fprintf([calling_function ' : ' 'Mixed \n'])
    t = lower_lim + indexMaxChange - 1;
end

end