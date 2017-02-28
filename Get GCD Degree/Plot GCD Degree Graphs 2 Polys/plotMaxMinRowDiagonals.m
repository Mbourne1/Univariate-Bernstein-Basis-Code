function plotMaxMinRowDiagonals(vMaxDiagR1,vMinDiagR1, limits)
%
% % Inputs
%
% vMaxDiagR1 :
%
% vMinDiagR1 : 
%
% limits :

lower_lim = limits(1);
upper_lim = limits(2);

x = lower_lim : 1 : upper_lim;

% Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
figure('name',figure_name)

vRatio_MaxMin_Diagonals_R = log10(vMinDiagR1)./log10(vMaxDiagR1);

plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
xlim([1 upper_lim]);
vline(lower_lim,'b','');
vline(upper_lim,'b','');
grid on
hold on
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off

end



