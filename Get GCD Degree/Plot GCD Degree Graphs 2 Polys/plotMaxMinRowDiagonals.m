function plotMaxMinRowDiagonals(vMaxDiagR1,vMinDiagR1, limits_k, limits_t)
%
% % Inputs
%
% vMaxDiagR1 : (Vector)
%
% vMinDiagR1 : (Vector) 
%
% limits_k : [Int Int]
%
% limits_t : [Int Int]

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

% 
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

% 
x = lowerLimit_k : 1 : upperLimit_k;

% Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
figure('name',figure_name)
vRatio_MaxMin_Diagonals_R = log10(vMinDiagR1)./log10(vMaxDiagR1);
plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
xlim([1 upperLimit_k]);
vline(lowerLimit_t,'b','');
vline(upperLimit_t,'b','');
grid on
hold on
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off

end



