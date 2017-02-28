function plotResiduals(vMinimumResiduals, my_limits, k_limits)
%
% % Inputs
%
% vMinimumResiduals : Vector of residuals
%
% limits : [lowerLimit upperLimit]


lowerLimit = my_limits(1);
upperLimit = my_limits(2);
vec_x = lowerLimit:1:upperLimit;

% % Plot Residuals
global SETTINGS
figure_name = sprintf('%s : Residuals of %s',mfilename, SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);

hold on

plot(vec_x,log10(vMinimumResiduals),'-s','DisplayName','Residuals by SVD')
ylabel('log_{10} Residuals')
xlabel('k')
legend(gca,'show');
vline(k_limits(1));
vline(k_limits(2));
hold off


end





