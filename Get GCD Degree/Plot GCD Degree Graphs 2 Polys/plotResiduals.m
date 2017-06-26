function plotResiduals(vMinimumResiduals, limits_k, limits_t, rank_range)
%
% % Inputs
%
% vMinimumResiduals : (Vector) Vector of residuals
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]
%
% rank_range


lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% % Plot Residuals
global SETTINGS
figure_name = sprintf('%s : Residuals of %s', mfilename, SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);

hold on
vec_x = lowerLimit_k:1:upperLimit_k;
plot(vec_x, log10(vMinimumResiduals), '-s', 'DisplayName', 'Residuals by SVD')
ylabel('log_{10} Residuals')
xlabel('k')
legend(gca,'show');

hline([rank_range mean(rank_range)],{'-r','-r','-b'})
vline(limits_t,{'r','r'})


hold off


end





