function plotRowDiagonals(arr_RowDiagonals, limits_k, limits_t)
%
% % Inputs
%
% arr_RowDiagonals : (Array of Vectors) containing the diagonals of the matrices R_{k}, 
% from the QR decomposition of S_{k} for k = lower_lim:upper_lim 
%
% limits_k : [Int Int] : Limits on the possible values of k
%
% limits_t : [Int Int] : Limits 


global SETTINGS

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(1);

% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
figure_name = sprintf('%s : Row Sum Norm',mfilename);
figure('name',figure_name)
hold on

for i = 1 : 1 : length(arr_RowDiagonals)
    
    % Get set of diagonals
    vec_RowDiags = arr_RowDiagonals{i};
    vec_i = i.*ones(length(vec_RowDiags));
    plot(vec_i, log10(vec_RowDiags) ,'*')

end

xlabel('k')
ylabel(sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_BUILD_METHOD))
title(sprintf('Diagonals of R1 from QR decomposition of %s',SETTINGS.SYLVESTER_BUILD_METHOD));

xlim([lowerLimit_k upperLimit_k]);
vline(lowerLimit_t, 'b', '');
vline(upperLimit_t, 'b', '');

hold off
end