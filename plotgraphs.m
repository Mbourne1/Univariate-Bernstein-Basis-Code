global SETTINGS
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        x = lower_lim:1:upper_lim;
        
        
        % Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
        figure('name',figure_name)
        vRatio_MaxMin_Diagonals_R = vMaxDiagR1./vMinDiagR1;
        plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
        hold on
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
        
        
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Norms',mfilename);
        figure('name',figure_name)
        vRatio_MaxMin_RowNorm_R = vMaxRowNormR1 ./ vMinRowNormR1;
        plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
        hold on
        legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
        hold off
        
        
        % Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
        figure_name = sprintf('%s : RowNorm',mfilename);
        figure('name',figure_name)
        plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
        xlabel('k')
        ylabel('Normalised Row Sums of R1 in S_{k}')
        title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
        
        
        % %
        figure_name = sprintf('%s : RowNorm',mfilename);
        figure('name',figure_name)
        plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
        xlabel('k')
        ylabel('Normalised Diagonals of R1 in S_{k}')
        title(['Normalised Diagonals in R1 matrix from the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
        
    case 'n'
        % Do Nothing
    otherwise
        error('err');
end