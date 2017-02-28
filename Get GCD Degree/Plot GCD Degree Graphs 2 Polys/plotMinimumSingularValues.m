
function plotMinimumSingularValues(vMinimumSingularValues, limits)
%
% % Inputs 
% 
% vMinimumSingularValues :
%
% 

lowerLimit = limits(1);
upperLimit = limits(2);

global SETTINGS

figure_name = sprintf('%s : Minimum Singular Values of %s',mfilename,SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name);
hold on
xlim([1 +inf]);
x_vec = lowerLimit:1:upperLimit;
plot(x_vec,log10(vMinimumSingularValues),'-s','DisplayName','Singular Values')
ylabel('log \sigma(k)')
xlabel('k')
legend(gca,'show');
hold off

end
