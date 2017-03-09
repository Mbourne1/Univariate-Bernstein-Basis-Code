function plotRowNorms(arr_RowNorms, myLimits, k_limits)
% 
% % Inputs
%
% arr_RowNorms 
%
% myLimits : These are defined by me
%
% k_limits : These are precomputed


global SETTINGS

% Set my limits
myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

% Set limits
lowerLimit = k_limits(1);
upperLimit = k_limits(2);

figure_name = sprintf('%s : Diag Norm',mfilename);
figure('name',figure_name)
hold on

for i = 1:1:length(arr_RowNorms)

    k = myLowerLimit + (i-1);
    
    % get vector of row norms
    vec_RowNorms = arr_RowNorms{i};
    vec_k = k.*ones(length(vec_RowNorms));
    
    plot(vec_k, log10(vec_RowNorms),'*');

    
end

vline(lowerLimit);
vline(upperLimit);

hold off

xlabel('k')
ylabel('log10 Row Norm of R1 from QR decomposition of S_{k}')
title(sprintf('log10 Row Norm of R1 from the QR decomposition of each subresultant %s', SETTINGS.SYLVESTER_BUILD_METHOD));
xlim([myLowerLimit, myUpperLimit]);

vline(myLowerLimit,'b','');
vline(myUpperLimit,'b','');
hold off
end