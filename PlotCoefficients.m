
function [] = PlotCoefficients(arrPolys, arrLabels, arrStyles)
% Plot coefficients of the set of polynomials stored in the array
%
% % Inputs
%
% arrPolys : (Array of vectors)
%
% arrLabels : (Array of Strings)

global SETTINGS
if SETTINGS.PLOT_GRAPHS
    
    nPolys = length(arrPolys);
    figure()
    hold on
    
    vDegree = zeros(nPolys,1);
    
    for i = 1 : 1 : nPolys
        
        mystyle = arrStyles{i};
        fx = arrPolys{i};
        
        name = arrLabels{i};
        
        vDegree(i) = GetDegree(fx);
        
        vec_x = 0 : 1 : vDegree(i);
        
        bool_log = true;
        if (bool_log == true)
            fx = log10(abs(fx));
        end
        
        plot(vec_x, fx, mystyle, 'DisplayName',name,'LineWidth',2)
        
    end
    
    xlim([1, max(vDegree)]);
    
    % Labels and legends
    xlabel('$i$ : Coefficient Index','Interpreter','latex', 'FontSize', 20);
    ylabel('$\log_{10} \left( \Re \right)$', 'Interpreter', 'latex', 'FontSize',20);
    l = legend(gca,'show');
    set(l,{'Interpreter','FontSize','Location'},{'latex',20, 'southwest'});
    hold off


    % Figure size and location
    myplot = gca;
    myval_side = 0.10;
    myval_base = 0.08;
    set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
    set(gcf, 'Position', [100, 100, 710, 650])
    
    box on
    grid on
    
end



end
