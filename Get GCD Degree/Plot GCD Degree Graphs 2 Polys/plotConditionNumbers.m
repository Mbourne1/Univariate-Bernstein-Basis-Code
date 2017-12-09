function [] = plotConditionNumbers(vConditionNumbers, limits_k, limits_t)
%
% % Inputs
%
% vConditionNumbers : (Vector of Floats)
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

global SETTINGS

figure_name = sprintf('%s : Condition Numbers of %s',mfilename,SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);
hold on
try
xlim([lowerLimit_k upperLimit_k]);
catch
end
x_vec = lowerLimit_k : 1 : upperLimit_k;


plot_name = strcat('$\kappa$ : BoolPreproc : ' , num2str(SETTINGS.BOOL_ALPHA_THETA));
plot(x_vec, log10(vConditionNumbers), '-s','DisplayName', plot_name)

vline(limits_t,{'r','r'})


% Plot Labels and legends
ylabel('$\log_{10} \left( \kappa_{k} \right)$', 'Interpreter', 'latex', 'FontSize',20)
xlabel('$k$', 'Interpreter', 'latex','FontSize',20)

l = legend(gca, 'show');
set(l,{'Interpreter', 'FontSize','Location'},{'latex', 20, 'southwest'});


% Appearance
grid on
box on


% Setting Size of plot within window
myplot = gca;
myval_side = 0.10;
myval_base = 0.08;
set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])


myWindow = gcf;
windowWidth = 700;
windowHeight = 610;
set(myWindow, 'Position', [ 100 100 windowWidth windowHeight])


hold off