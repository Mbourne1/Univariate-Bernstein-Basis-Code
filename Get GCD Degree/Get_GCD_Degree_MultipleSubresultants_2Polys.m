function [t] = Get_GCD_Degree_MultipleSubresultants_2Polys(vMetric, deg_limits)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.
%
% % Inputs
%
% vMetric : vector of values used to determine degree of GCD. These values
% may be 
%       'Minimum Singular values of S_{k}'
%       'Min/Max Row Diagonals of R_{k}'
%       'Min/Max Row Norms of R_{k}'
%
% deg_limits : Limits of degree of GCD


% Intialise global settings
global SETTINGS

% Get the function which called this function.
[St,~] = dbstack();
calling_function = St(2).name;

lowerLimit = deg_limits(1);
upperLimit = deg_limits(2);


% Get the maximum change in singular values and the index at which the
% maximum change occured.
[maxChangeMetric, indexMaxChange] = Analysis(vMetric);

display([mfilename ' : ' calling_function ' : ' sprintf('Max Change : %2.4f ', maxChangeMetric)]);
display([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %2.4f ',SETTINGS.THRESHOLD)]);


if  abs(maxChangeMetric) < SETTINGS.THRESHOLD
    
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMetricValue = log10(mean(vMetric));
    
    % %
    % %
    % %
    if(SETTINGS.PLOT_GRAPHS)
        
        figure_name = sprintf([mfilename ' : ' calling_function ': Singular Values of %s'],SETTINGS.SYLVESTER_BUILD_METHOD);
        figure('name',figure_name)
        plot(log10(vMetric),'DisplayName','Singular Values');
        hold on
        mu = avgMetricValue;
        hline = refline([0 mu]);
        hline.Color = 'r';
        vline(lowerLimit,'b','');
        vline(lowerLimit,'b','');
        hold off
        
        
    end
    
    
    if  avgMetricValue < SETTINGS.THRESHOLD_RANK
        % If all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        t = upperLimit;
        
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
    t = lowerLimit + indexMaxChange - 1;
    
    fprintf([mfilename ' : ' calling_function ' : ' 'Mixed \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
    
    
end

end