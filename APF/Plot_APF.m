
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Residuals',mfilename);
        figure('name',figure_name)
        title('plotting residual')
        hold on
        plot(1:1:length(residual),(residual));
        
    case 'n'
    otherwise
        error('err')
end
