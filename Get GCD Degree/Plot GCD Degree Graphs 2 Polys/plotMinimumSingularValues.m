
function plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t)
%
% % Inputs 
% 
% vMinimumSingularValues : (Vector)
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

global SETTINGS

figure_name = sprintf('%s : Minimum Singular Values of %s',mfilename,SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);
hold on
xlim([1 +inf]);

x_vec = lowerLimit_k:1:upperLimit_k;

plot(x_vec,log10(vMinimumSingularValues),'-s','DisplayName','Singular Values')
ylabel('log \sigma(k)')
xlabel('k')
legend(gca,'show');

vline(lowerLimit_t);
vline(upperLimit_t);

hold off

end
