function plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits)
% 
% % Inputs
%
% vMaxRowNormR1
%
% vMinRowNormR1
%
% limits

global SETTINGS

[lower_lim] = limits(1);
[upper_lim] = limits(2);

x = lower_lim:1:upper_lim;

% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Norms',mfilename);
figure('name',figure_name)

vRatio_MaxMin_RowNorm_R = vMaxRowNormR1 ./ vMinRowNormR1;

plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
hold on

legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
title(sprintf('Max:Min Row Norms of Rows in R1 from the QR Decomposition of %s', SETTINGS.SYLVESTER_BUILD_METHOD));
xlim([1 upper_lim]);

vline(lower_lim,'b','');
vline(upper_lim,'b','');

hold off
end
