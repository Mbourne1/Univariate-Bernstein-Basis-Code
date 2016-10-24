% 1. Row Diagonals of matrix R_{k} from QR decomposition of S_{k}
% 2. Row Norms of matrix R_{k} from QR decompositon of S_{k}
%
% 4. Diagonals of R1


global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        x = lower_lim_comp:1:upper_lim_comp;
        
        % 1.
         
        % Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
        figure('name',figure_name)
        vRatio_MaxMin_Diagonals_R = vMaxDiagR1./vMinDiagR1;
        
        plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
        xlim([1 upper_lim_comp]);
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold on
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
        
        % 2.
        
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Norms',mfilename);
        figure('name',figure_name)
        vRatio_MaxMin_RowNorm_R = vMaxRowNormR1 ./ vMinRowNormR1;
        plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
        hold on
        legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
         xlim([1 upper_lim_comp]);
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        
        % 3.
        
        % Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
        figure_name = sprintf('%s : RowNorm',mfilename);
        figure('name',figure_name)
        plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
        xlabel('k')
        ylabel('Normalised Row Sums of R1 in S_{k}')
        title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
         xlim([1 upper_lim_comp]);
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        
        % 4.
        % %
        figure_name = sprintf('%s : RowNorm',mfilename);
        figure('name',figure_name)
        plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
        xlabel('k')
        ylabel('Normalised Diagonals of R1 in S_{k}')
        title(['Normalised Diagonals in R1 matrix from the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        xlim([1 upper_lim_comp]);
        vline(lower_lim,'b','');
        vline(upper_lim,'b','');
        hold off
        
        %5.
        % % Plot Residuals
        figure_name = sprintf('%s : Residuals',mfilename);
        figure('name',figure_name);
        hold on
        plot(log10(vMinimumResidual_SVD),'-s','DisplayName','SVD')
        plot(log10(vMinimumResidual_QR),'-s','DisplayName','QR')
        ylabel('log r(k)')
        xlabel('k')
        legend(gca,'show');
        hold off
        
        
    case 'n'
        % Do Nothing
    otherwise
        error('err');
end