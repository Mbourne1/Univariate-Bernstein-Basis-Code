function t = GetRankOneSubresultant(DTQ)
% Given the vector of values from either minimum singular values or max:min
% R diagonals.
% Get the rank, where only one subresultant exists.


global SETTINGS

[St,~] = dbstack();
calling_function = St(2).name;

vMinSingVal = svd(DTQ);

% Only one subresultant
fprintf('Only one subresultant exists. \n')
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

[deltaSingularValues,~] = Analysis(vMinSingVal);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.
if deltaSingularValues < SETTINGS.THRESHOLD
    
    % The subresultant is of full rank, in which case t = 0
    t = 0;
    fprintf([calling_function ' : ' 'The only Subresultant S_{1} appears to be of NonSingular. \n']);
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = 1;
    fprintf([calling_function ' : ' 'The only Subresultant S_{1} appears to be Singular \n']);
    return
    
end
end