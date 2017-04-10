function plotResiduals(vMinimumResiduals, limits_k, limits_t)
%
% % Inputs
%
% vMinimumResiduals : (Vector) Vector of residuals
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);
vec_x = lowerLimit_k:1:upperLimit_k;

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);


% % Plot Residuals
global SETTINGS
figure_name = sprintf('%s : Residuals of %s', mfilename, SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);

hold on

plot(vec_x, log10(vMinimumResiduals), '-s', 'DisplayName', 'Residuals by SVD')
ylabel('log_{10} Residuals')
xlabel('k')
legend(gca,'show');
vline(lowerLimit_t);
vline(upperLimit_t);
hold off


end





