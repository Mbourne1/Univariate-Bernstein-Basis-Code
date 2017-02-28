
function plotRowDiagonals(arr_RowDiagonals, k_limits)
%
% % Inputs
%
% arr_RowDiagonals : Array containing the diagonals of the matrices R_{k}, 
% from the QR decomposition of S_{k} for k = lower_lim:upper_lim 
%
% k_limits : [lower_lim upper_lim] : Limits on the possible values of k

global SETTINGS

lower_lim = k_limits(1);
upper_lim = k_limits(2);

% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
figure_name = sprintf('%s : Row Sum Norm',mfilename);
figure('name',figure_name)
hold on

for i = 1:1:length(arr_RowDiagonals)
    
    % Get set of diagonals
    vec_RowDiags = arr_RowDiagonals{i};
    vec_i = i.*ones(length(vec_RowDiags));
    
    plot(vec_i, log10(vec_RowDiags) ,'*')
end

xlabel('k')

ylabel(sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_BUILD_METHOD))
title(sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_BUILD_METHOD));

xlim([1 upper_lim]);
vline(lower_lim,'b','');
vline(upper_lim,'b','');

hold off
end