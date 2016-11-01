function t = Get_GCD_Degree_OneSubresultant(vSingularValues)
% Given the vector of values from either minimum singular values or max:min
% R diagonals.
% Get the rank, where only one subresultant exists.
%
% Inputs 
%
% vSingularValues : Vector of Singular Values.
%
%
% Outputs.
%
% t : Computed degree of the GCD


global SETTINGS

% Get the calling function
[St,~] = dbstack();
calling_function = St(2).name;



% Only one subresultant
fprintf([calling_function ' : ' mfilename ' : ' 'Only one subresultant exists. \n'])

% Plot the singular values
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf([calling_function ' : Singular values of S_{1}']);
        figure('name',figure_name)
        hold on
        title('Singular values of S_{1}')
        plot(log10(vMinSingVal))
        hold off
    case 'n'
end

[deltaSingularValues,~] = Analysis(vSingularValues);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.
if deltaSingularValues < SETTINGS.THRESHOLD
    
    % The subresultant is of full rank, in which case t = 0
    t = 0;
    fprintf([mfilename ' : ' calling_function ' : ' 'The only Subresultant S_{1} appears to be of full rank. \n']);
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = 1;
    fprintf([mfilename ' : ' calling_function ' : ' 'The only Subresultant S_{1} appears to be rank deficient \n']);
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
    return
    
end
end