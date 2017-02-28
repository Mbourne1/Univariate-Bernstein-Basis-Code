
function plotSingularValues(arr_SingularValues, limits)
%
% % Inputs
%
% arr_SingularValues
%
% limits

lowerLimit = limits(1);
upperLimit = limits(2);
nSubresultants = upperLimit - lowerLimit + 1;

global SETTINGS

figure_name = sprintf('%s : All Singular Values of %s ',mfilename, SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name', figure_name);

hold on
for i = 1: 1: nSubresultants
   
    k = lowerLimit + (i-1);
    
    vSingularValues = arr_SingularValues{i};
    
    vec_k = k.*ones(length(vSingularValues));
    
    plot(vec_k, log10(vSingularValues), '*')
    
end
hold off


end