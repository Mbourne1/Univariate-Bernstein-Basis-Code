function plotResiduals(vMinimumResiduals, myLimits, k_limits)
%
% % Inputs
%
% vMinimumResiduals : Vector of residuals
%
% limits : [lowerLimit upperLimit]


myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);
vec_x = myLowerLimit:1:myUpperLimit;

lowerLimit = k_limits(1);
upperLimit = k_limits(2);


% % Plot Residuals
global SETTINGS
figure_name = sprintf('%s : Residuals of %s', mfilename, SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);

hold on

plot(vec_x, log10(vMinimumResiduals), '-s', 'DisplayName', 'Residuals by SVD')
ylabel('log_{10} Residuals')
xlabel('k')
legend(gca,'show');
vline(lowerLimit);
vline(upperLimit);
hold off


end





