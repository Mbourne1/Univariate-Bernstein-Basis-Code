function [t, alpha, theta, GM_fx, GM_gx, rank_range] = ...
    Get_GCD_Degree_2Polys(fx, gx, limits_t, rank_range)
% GetGCD_Degree_2Polys(fx, gx, t_limits)
%
% Get degree t of the AGCD d(x) of input polynomials f(x) and g(x)
%
%
% % Inputs.
%
% fx : (Vector) coefficients of polynomial f, expressed in Bernstein Basis
%
% gx : (Vector) coefficients of polynomail g, expressed in Bernstein Basis
%
% t_limits : [(Int) (Int)] Set the upper and lower bound of the degree of the
% GCD of polynomials f(x) and g(x). Usually used when using o_roots() where
% prior information is known about the upper and lower bound due to number
% of distinct roots.
%
% rank_range : [Float Float]
%
% % Outputs.
%
% t : (Int) Degree of GCD of f(x) and g(x)
%
% alpha : (Float) Optimal value of alpha
%
% theta : (Float) Optimal value of theta
%
% deg_limits : [(Int) (Int)]
%
% rank_range : [Float Float]

global SETTINGS

% Get degree of polynomail f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Set my limits to either be equal to the precalculated limits, or an
% extended range.
limits_k = [1 min(m,n)];

% If the number of distinct roots in f(x) is one, then the degree of the
% GCD of f(x) and f'(x) = m-1 = n.
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

%
rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

% Get the number of subresultants which must be constructed.
nSubresultants = upperLimit_k - lowerLimit_k +1 ;

% %
% Initialisation of some vectors

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(nSubresultants, 1);
vTheta    =   zeros(nSubresultants, 1);
vGM_fx    =   zeros(nSubresultants, 1);
vGM_gx    =   zeros(nSubresultants, 1);
vConditionSk = zeros(nSubresultants, 1);
vConditionCf = zeros(nSubresultants, 1);
vConditionCg = zeros(nSubresultants, 1);

% Initialise some cell arrays
arr_SingularValues = cell(nSubresultants, 1);
arr_Sk = cell(nSubresultants, 1);
arr_R  = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);



% For each subresultant $S_{k} k = lowerLimt ... upperLimit$
for k = lowerLimit_k : 1 : upperLimit_k
    
    % Set index i
    i = k - lowerLimit_k + 1;
    
    % Get Geometric means, \alpha, and \theta
    [vGM_fx(i), vGM_gx(i), vAlpha(i), vTheta(i)] = Preprocess_2Polys(fx, gx, k);
    
    % 18/04/2016
    % Given the previous geometric mean of f(x) calculate the new geometric
    % mean by my new method. 
    if (i > 1)
        GM_fx_test = GetGeometricMeanFromPrevious(fx , vGM_fx(i-1) , m , n-k);
        GM_gx_test = GetGeometricMeanFromPrevious(gx , vGM_gx(i-1) , n , m-k);
    end
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx./ vGM_fx(i);
    gx_n = gx./ vGM_gx(i);
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = GetWithThetas(fx_n, vTheta(i));
    gw = GetWithThetas(gx_n, vTheta(i));
    
    plotCoefficients(fx, fw);
    plotCoefficients(gx, vAlpha(i).*gw);
    
    % Build the k-th subresultant
    arr_Sk{i} = BuildSubresultant_2Polys(fw, vAlpha(i).*gw, k);
    
    % Get condition of the sylvester subresultant matrix
    vConditionSk(i) = cond(arr_Sk{i});
    
    % Get condition of first partiton
    % vConditionCf(i) = cond(arr_Sk{i}(:, 1:n-k+1));
    % vConditionCg(i) = cond(arr_Sk{i}(:, n-k+2:end));
    
    % Get condition of second partition
    
    
    % Get singular values of S_{k}
    arr_SingularValues{i} = svd(arr_Sk{i});
    
    % Get the QR decomposition
    [~, arr_R{i}] = qr(arr_Sk{i});
    [~,c] = size(arr_R{i});
    arr_R1{i} = arr_R{i}(1:c,1:c);
    
    
    
end

plotConditionNumbers(vConditionSk, limits_k, limits_t)
%plotConditionNumbers(vConditionCf, limits_k, limits_t)
%plotConditionNumbers(vConditionCg, limits_k, limits_t)


fprintf(sprintf('Metric used to compute degree of GCD : %s \n', SETTINGS.RANK_REVEALING_METRIC));

% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

switch SETTINGS.RANK_REVEALING_METRIC
    
    case 'R1 Row Norms'
        
        % Preallocate vectors
        vMaxRowNormR1 = zeros(nSubresultants, 1);
        vMinRowNormR1 = zeros(nSubresultants, 1);
        arr_R1_RowNorms = cell(nSubresultants, 1);
        
        % Get maximum and minimum row norms of rows of R1.
        for i = 1:1: nSubresultants
            
            arr_R1_RowNorms{i} = sqrt(sum(arr_R1{i}.^2,2)) ./ norm(arr_R1{i});
            
            
            vMaxRowNormR1(i) = max(arr_R1_RowNorms{i});
            vMinRowNormR1(i) = min(arr_R1_RowNorms{i});
        end
        
        % Plot graphs
        if (SETTINGS.PLOT_GRAPHS)
            plotRowNorms(arr_R1_RowNorms, limits_k, limits_t);
            plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
        end
        
        vMetric = log10(vMinRowNormR1./vMaxRowNormR1);
        
    case 'R1 Row Diagonals'
        
        % Initialise vectors to store max and min diagonal entries of R1
        vMaxDiagR1 = zeros(nSubresultants, 1);
        vMinDiagR1 = zeros(nSubresultants, 1);
        arr_DiagsR1 = cell(nSubresultants, 1);
        
        % Get maximum and minimum row diagonals of R1
        for i = 1:1:nSubresultants
            arr_DiagsR1{i} = diag(abs(arr_R1{i}));
            
            vMaxDiagR1(i) = max(arr_DiagsR1{i});
            vMinDiagR1(i) = min(arr_DiagsR1{i});
            
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot Graphs
            plotRowDiagonals(arr_DiagsR1, limits_k, limits_t);
            plotMaxMinRowDiagonals(vMaxDiagR1,vMinDiagR1, limits_k, limits_t, rank_range);
        end
        
        vMetric = log10(vMinDiagR1./vMaxDiagR1);
        
    case 'Minimum Singular Values'
        
        % Get the minimal singular value from S_{k}
        
        % Initialise a vector to store minimimum singular value of each
        % S_{k}
        vMinimumSingularValues = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            
            % Get minimum Singular value of S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
        end
        
        vMetric = log10(vMinimumSingularValues);
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot Graphs
            plotSingularValues(arr_SingularValues, limits_k, limits_t);
            plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t, rank_range)
        end
    
    case 'Max/Min Singular Values'
        
        % Get the minimal singular value from S_{k}
        
        % Initialise a vector to store minimimum singular value of each
        % S_{k}
        vMinimumSingularValues = zeros(nSubresultants, 1);
        vMaximumSingularValues = zeros(nSubresultants, 1);
        
        for i = 1:1:nSubresultants
            
            % Get minimum Singular value of S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
            vMaximumSingularValues(i) = max(arr_SingularValues{i});
        end
        
        vMetric = log10(vMinimumSingularValues) - log10(vMaximumSingularValues);
        
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot Graphs
            plotSingularValues(arr_SingularValues, limits_k, limits_t);
            plotMaxMinSingularValues(vMaximumSingularValues, vMinimumSingularValues, limits_k, limits_t, rank_range)
        end
        
        
    case 'Residuals'
        
        % Initialise vectors to store residuals
        vMinimumResidual_QR  = zeros(nSubresultants,1);
        vMinimumResidual_SVD = zeros(nSubresultants,1);
        
        % Get the minimal residuals for each subresultant S_{k} by QR and
        % SVD
        
        for i = 1:1:nSubresultants
            
            vMinimumResidual_QR(i) = GetMinimalDistance(arr_Sk{i},'QR');
            vMinimumResidual_SVD(i) = GetMinimalDistance(arr_Sk{i},'SVD');
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot Graphs
            plotResiduals(vMinimumResidual_QR, limits_k, limits_t, rank_range);
            %plotResiduals(vMinimumResidual_SVD, k_limits);
        end
        
        vMetric = log10(vMinimumResidual_QR);
        
    otherwise
        error('%s is not a valid rank revealing metric', SETTINGS.RANK_REVEALING_METRIC)
        
end


LineBreakLarge();

% If only one subresultant exists, use an alternative method.

if (upperLimit_t == lowerLimit_t)
        
    t = upperLimit_t;
    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);

    
elseif (upperLimit_k == lowerLimit_k )
    
    % Get the metric to compute the degree of the GCD using only one
    % subresultant matrix.
    vMetric = log10(arr_SingularValues{1});
    
    % Compute the degree of the GCD.
    t = Get_GCD_Degree_OneSubresultant_2Polys(vMetric, rank_range);
    

    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    
    
    
else
    
    [t] = Get_GCD_Degree_MultipleSubresultants_2Polys(vMetric, limits_k, limits_t, rank_range);
    
    rank_range_low = vMetric(t - (lowerLimit_k-1));
    i = t - (lowerLimit_k - 1);
    
    if i ~= (upperLimit_k )
        rank_range_high = vMetric(t - (lowerLimit_k-1) + 1);
    end
    
    rank_range = [rank_range_low rank_range_high];
    
    % Output all subresultants, all optimal alphas, all optimal thetas and all
    % geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
    
    
    alpha = vAlpha(i);
    theta = vTheta(i);
    GM_fx = vGM_fx(i);
    GM_gx = vGM_gx(i);
    
end




end


function [] = plotCoefficients(fx,fw)

m = GetDegree(fx);

figure()
hold on
vec_x = 0:1:m;
plot(vec_x, log10(fx),'r-s','DisplayName','f(x)');
plot(vec_x, log10(fw),'b-o','DisplayName','f(\omega)');
legend(gca,'show');
xlabel('i')
ylabel('log_{10} (\Re)')
hold off

end







