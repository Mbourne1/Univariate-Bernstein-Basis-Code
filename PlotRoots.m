global PLOT_GRAPHS

% given the roots of f and g, plot them on a line
switch PLOT_GRAPHS
    case 'n'
        % Dont plot graphs
    case 'y'
        figure('name','Exact roots of f(x), g(x) and d(x)')
        hold on
        title('Roots of f and g on the real interval')
        scatter(f_roots(:,1),ones(size(f_roots(:,1))),'s','DisplayName','Roots of f(x)')
        try
            scatter(g_roots(:,1),ones(size(g_roots(:,1))),'x','DisplayName','Roots of g(x)')
        catch
            fprintf('could not plot exact roots of g\n')
        end
        try
            scatter(d_roots(:,1),ones(size(d_roots(:,1))),'o','DisplayName','Roots of d(x)')
        catch
            fprintf('Could not plot exact roots of d.\n')
        end
        xlabel('Real')
        legend(gca,'show')
        hold off
    otherwise
        error('error PLOT_GRAPH is either y or n')
end