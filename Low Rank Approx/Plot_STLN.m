
if(SETTINGS.PLOT_GRAPHS)
        figure_name = sprintf([mfilename ' : ' 'Residuals']);
        figure('name',figure_name)
        hold on
        xlim([1 +inf]);
        plot(log10(condition),'-s');
        hold off
        
end