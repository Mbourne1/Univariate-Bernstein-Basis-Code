function plotMaxMinRowDiagonals(vMaxDiagR1, vMinDiagR1, limits_k, limits_t, rank_range)
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
% rank_range : [Float Float]

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


% 
x_vec = lowerLimit_k : 1 : upperLimit_k;

% Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
figure('name',figure_name)

vRatio_MaxMin_Diagonals_R = vMinDiagR1./vMaxDiagR1;

plot(x_vec, log10(vRatio_MaxMin_Diagonals_R),'red-s');
xlim([lowerLimit_k upperLimit_k]);

hline([rank_range mean(rank_range)],{'-r','-r','-b'})

vline(limits_t,{'-r','-r'})

grid on
hold on
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off

end



