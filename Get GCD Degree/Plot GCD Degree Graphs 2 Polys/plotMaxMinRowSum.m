function plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
% 
% % Inputs
%
% vMaxRowNormR1 : (Vector)
%
% vMinRowNormR1 : (Vector)
%
% limits_k : [(Int) (Int)]
%
% limits_t : [(Int) (Int)]
%
% rank_range : [(Float) (Float)]
% 
% % Outputs



global SETTINGS

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
x_vec = lowerLimit_k : 1 : upperLimit_k;

% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Norms',mfilename);
figure('name',figure_name)

vRatio_MaxMin_RowNorm_R = vMaxRowNormR1 ./ vMinRowNormR1;
plot(x_vec, log10(vRatio_MaxMin_RowNorm_R),'red-s');
hold on

legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
title(sprintf('Max:Min Row Norms of Rows in R1 from the QR Decomposition of %s', SETTINGS.SYLVESTER_BUILD_METHOD));
%xlim([1 upperLimit_k]);

hline(rank_range, {'r','r'});
vline(limits_t,{'b','b'});


hold off
end
