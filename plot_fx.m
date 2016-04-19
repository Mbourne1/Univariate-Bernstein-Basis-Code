function [] = Plot_fx(fx,a,b,figure_name)
% Plot the Bernstein polynomial f(x) and its control points defined over
% the interval [a,b]


global PLOT_GRAPHS
switch PLOT_GRAPHS
    case 'y'
        
        % Get the control points of f(x)
        Pk = GetControlPoints(a,b,fx);
        
        % Initialise the column vector of x values
        x = linspace(a,b,100)';
        
        % Get number of entries in x vector
        nEntries_x = size(x,1);
        
        % Initialise a y vector of the same size as x
        y = zeros(nEntries_x,1);
        
        % Set the y ordinate values.
        for i = 1:1:nEntries_x
            
            % Evaluate f(x) at x_{i}
            y(i) = Bernstein_Evaluate(fx,x(i));
            
        end
        
        
        figure('name',figure_name)
        plot(x,y,'-');
        hold on
        scatter(Pk(:,1),Pk(:,2));
        hold off
    case 'n'
    otherwise
        error('err');
end
end