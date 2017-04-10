function plotRowNorms(arr_RowNorms, limits_k, limits_t)
% 
% % Inputs
%
% arr_RowNorms : (Array of Vectors) 
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]

global SETTINGS

% Set my limits
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% Set limits
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

figure_name = sprintf('%s : Diag Norm',mfilename);
figure('name',figure_name)
hold on

for i = 1:1:length(arr_RowNorms)

    k = lowerLimit_k + (i-1);
    
    % get vector of row norms
    vec_RowNorms = arr_RowNorms{i};
    vec_k = k.*ones(length(vec_RowNorms));
    
    plot(vec_k, log10(vec_RowNorms),'*');

    
end

xlabel('k')
ylabel('log10 Row Norm of R1 from QR decomposition of S_{k}')
title(sprintf('log10 Row Norm of R1 from the QR decomposition of each subresultant %s', SETTINGS.SYLVESTER_BUILD_METHOD));
xlim([lowerLimit_k, upperLimit_k]);

vline(lowerLimit_t, 'b', '');
vline(upperLimit_t, 'b', '');

hold off
end