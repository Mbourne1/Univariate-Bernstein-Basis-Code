
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        figure_name = sprintf([mfilename ' : ' 'Residuals']);
        figure('name',figure_name)
        hold on
        xlim([1 +inf]);
        plot(log10(condition),'-s');
        hold off
        
    case 'n'
        
    otherwise
        
        error('SETTINGS.PLOT_GRAPHS must be either y or n')
        
end