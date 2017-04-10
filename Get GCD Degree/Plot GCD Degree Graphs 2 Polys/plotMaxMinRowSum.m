function plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t)
% 
% % Inputs
%
% vMaxRowNormR1 :
%
% vMinRowNormR1 :
%
% limits_k :
%
% limits_t :
%
% % Outputs
%


global SETTINGS

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

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
xlim([1 upperLimit_k]);

vline(lowerLimit_t,'b','');
vline(upperLimit_t,'b','');

hold off
end
