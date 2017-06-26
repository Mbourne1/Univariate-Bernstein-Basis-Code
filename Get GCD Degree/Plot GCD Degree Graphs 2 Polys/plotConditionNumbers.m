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
plot(x_vec, log10(vConditionNumbers), '-s','DisplayName','Condition Numbers')

vline(limits_t,{'r','r'})


% Plot Labels
ylabel('log (\kappa{k})')
xlabel('k')
legend(gca,'show');



hold off