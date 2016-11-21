
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Condition APF',mfilename);
        figure('name',figure_name)
        title('Plotting condition')
        hold on
        plot(1:1:length(condition),(condition));
        
    case 'n'
    otherwise
        error('err')
end
