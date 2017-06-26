function [] = MaxMinEntries(m, n_k)
% MaxMinEntriesDT(m,n_k)
%
% Each coefficient a_{i} appears in n-k+1 columns of the Sylvester matrix,
% and has two binomial coefficients in D^{-1}T(f,g).
% This experiment looks at the scaling effect of the two binomial
% coefficients of each a_{i} in each column j = 0,...,n-k.
%
% Inputs
%
% m : Degree of polynomial f(x)
%
% n : Degree of polynomial v(x)
%
% >> MaxMinEntriesDT(m,n_k)


% Get an array of Sylvester matrix formats
arrSylvesterMatrixType = {'DTQ', 'TQ', 'DT', 'T', 'DTQ Denominator Removed'};

% Get number of formats in the array
nFormats = length(arrSylvesterMatrixType);

% for each coefficient a_{i} i=0,...,m
arrScaling = cell(m+1,1);

% For each Sylvester subresultant format
for j = 1 : 1 : nFormats
    
    % For each of the coefficients a_{i}
    for i = 0 : 1 : m
        
        % Get format string
        subresultantFormat = arrSylvesterMatrixType{j};
        
        % Get vector of scaling of the coefficient a_{i} in each column of
        % the first partition of the subresultant matrix
        vScaling = GetScalingVector(m, n_k, i, subresultantFormat);
        
        % Add to array
        arrScaling{i+1,1} = vScaling;
        
    end
    
    
    
    % Begin the plotting
    figure_name = sprintf('%s : Scaling effect in %s',mfilename, subresultantFormat );
    figure('name',figure_name)
    hold on
    
    for i = 0 : 1 : m
        % Plot values for coefficient
        
        % Set name to appear in legend
        legend_name = sprintf('Coefficient : a_{%i}',i);
        
        plot(log10(arrScaling{i+1,1}),'DisplayName',legend_name)
        
    end
    
    legend(gca, 'show');
    ylabel('log_{10}(\Re)')
    xlabel('Column index')
    hold off
    
    
    
end



end



function vScaling = GetScalingVector(m, n_k, i, sylvesterMatrixFormat)
%
% % Inputs
%
% m : (Int)
%
% n_k : (Int)
%
% i : (Int)
%
% sylvesterMatrixFormat : (String)


vScaling = zeros(n_k + 1, 1);

switch sylvesterMatrixFormat
    case 'T'
        
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m,i);
        end
        
    case 'DT'
        
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m,i) ./ nchoosek(m + n_k, i + j);
        end
        
    case 'TQ'
        
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m, i) .* nchoosek(n_k, j);
        end
        
    case 'DTQ'
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m, i) .* nchoosek(n_k, j) ./ nchoosek(m+n_k, i + j);
        end
        
    case 'DTQ Denominator Removed'
        for j = 0 : 1 : n_k
            vScaling(j + 1) = nchoosek(m, i) .* nchoosek(n_k, j) ...
                ./ nchoosek(m+n_k, i + j) .* nchoosek(m+n_k, n_k);
        end
        
    otherwise
        error('err')
end
end