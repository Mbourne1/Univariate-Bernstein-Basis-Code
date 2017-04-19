function [t,alpha, theta,GM_fx,GM_gx] = ...
    Get_GCD_Degree_2Polys(fx, gx, t_limits)
% GetGCD_Degree_2Polys(fx, gx, k_limits)
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
% k_limits : [(Int) (Int)] Set the upper and lower bound of the degree of the
% GCD of polynomials f(x) and g(x). Usually used when using o_roots() where
% prior information is known about the upper and lower bound due to number
% of distinct roots.
%
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

global SETTINGS

% Get degree of polynomail f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Set my limits to either be equal to the precalculated limits, or an
% extended range.
k_limits = [1 min(m,n)];

% If the number of distinct roots in f(x) is one, then the degree of the
% GCD of f(x) and f'(x) = m-1 = n.
lowerLimit_k = k_limits(1);
upperLimit_k = k_limits(2);



% Get the number of subresultants which must be constructed.
nSubresultants = upperLimit_k - lowerLimit_k +1 ;

% %
% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(nSubresultants, 1);
vTheta    =   zeros(nSubresultants, 1);
vGM_fx    =   zeros(nSubresultants, 1);
vGM_gx    =   zeros(nSubresultants, 1);


% Initialise some cell arrays
arr_SingularValues = cell(nSubresultants, 1);
arr_Sk = cell(nSubresultants, 1);
arr_R  = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);

% For each subresultant $S_{k} k = lowerLimt ... upperLimit$
for k = lowerLimit_k : 1 : upperLimit_k
    
    i = k - lowerLimit_k + 1;
    
    [vGM_fx(i), vGM_gx(i), vAlpha(i), vTheta(i)] = Preprocess_2Polys(fx, gx, k);
    
    % 18/04/2016
    % Given the previous geometric mean of f(x) calculate the new geometric
    % mean by my new method
    
    if (i>1)
        GM_fx_test = GetGeometricMeanFromPrevious(fx , vGM_fx(i-1) , m , n-k);
        GM_gx_test = GetGeometricMeanFromPrevious(gx , vGM_gx(i-1) , n , m-k);
    end
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx./ vGM_fx(i);
    gx_n = gx./ vGM_gx(i);
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = GetWithThetas(fx_n, vTheta(i));
    gw = GetWithThetas(gx_n, vTheta(i));
    
    % Build the k-th subresultant
    arr_Sk{i} = BuildSubresultant_2Polys(fw, vAlpha(i).*gw, k);
    
    % Get singular values of S_{k}
    arr_SingularValues{i} = svd(arr_Sk{i});
    
    % Get the QR decomposition
    [~, arr_R{i}] = qr(arr_Sk{i});
    [r,c] = size(arr_R{i});
    arr_R1{i} = arr_R{i}(1:c,1:c);
    
    
    
end


if(SETTINGS.PLOT_GRAPHS)
    
    plot_extra_graphs = false;
    if (plot_extra_graphs)
        x_vec = lowerLimit_k:1:upperLimit_k;
        
        figure_name = sprintf('Geometric Mean of f(x) %s', SETTINGS.SYLVESTER_BUILD_METHOD);
        figure('name',figure_name)
        hold on
        plot(x_vec,log10(vGM_fx), '-s')
        xlabel('k')
        ylabel('log_{10} lamda_{k}')
        %xlabel('$\alpha$','Interpreter','LaTex')
        title(figure_name)
        hold off
        
        figure_name = sprintf('Geometric Mean of g(x) %s', SETTINGS.SYLVESTER_BUILD_METHOD);
        figure('name',figure_name)
        hold on
        plot(x_vec,log10(vGM_gx), '-s')
        xlabel('k')
        ylabel('log_{10} mu')
        title(figure_name)
        hold off
        
        figure_name = sprintf('Theta values in %s',SETTINGS.SYLVESTER_BUILD_METHOD);
        figure('name',figure_name)
        hold on
        plot(x_vec,log10(vTheta), '-s');
        xlabel('k')
        ylabel('log_{10} theta')
        title(figure_name)
        hold off
        
        figure_name = sprintf('Alpha values in %s', SETTINGS.SYLVESTER_BUILD_METHOD);
        title(figure_name)
        figure('name',figure_name)
        hold on
        plot(x_vec,log10(vAlpha), '-s')
        xlabel('k')
        ylabel('log_{10} alpha')
        title(figure_name)
        hold off
    end
end



fprintf(sprintf('Metric used to compute degree of GCD : %s \n',SETTINGS.RANK_REVEALING_METRIC));

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
            
            arr_R1_RowNorms{k} = sqrt(sum(arr_R1{i}.^2,2))./norm(arr_R1{i});
            
            vMaxRowNormR1(i) = max(arr_R1_RowNorms{i});
            vMinRowNormR1(i) = min(arr_R1_RowNorms{i});
        end
        
        % Plot graphs
        if (SETTINGS.PLOT_GRAPHS)
            plotRowNorms(arr_R1_RowNorms, k_limits, t_limits);
            plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, k_limits, t_limits)
        end
        
        metric = vMinRowNormR1./vMaxRowNormR1;
        
    case 'R1 Row Diagonals'
        
        % Initialise vectors to store max and min diagonal entries of R1
        vMaxDiagR1 = zeros(nSubresultants, 1);
        vMinDiagR1 = zeros(nSubresultants, 1);
        arr_DiagsR1 = cell(nSubresultants, 1);
        
        % Get maximum and minimum row diagonals of R1
        for i = 1:1:nSubresultants
            arr_DiagsR1{k} = diag(abs(arr_R1{i}));
            
            vMaxDiagR1(i) = max(arr_DiagsR1{i});
            vMinDiagR1(i) = min(arr_DiagsR1{i});
            
        end
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot Graphs
            plotRowDiagonals(arr_DiagsR1, k_limits, t_limits);
            plotMaxMinRowDiagonals(vMaxDiagR1,vMinDiagR1, k_limits, t_limits);
        end
        
        metric = vMinDiagR1./vMaxDiagR1;
        
    case 'Singular Values'
        
        % Get the minimal singular value from S_{k}
        
        % Initialise a vector to store minimimum singular value of each
        % S_{k}
        vMinimumSingularValues = zeros(nSubresultants,1);
        
        for i = 1:1:nSubresultants
            
            % Get minimum Singular value of S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
        end
        
        metric = vMinimumSingularValues;
        
        if (SETTINGS.PLOT_GRAPHS)
            % Plot Graphs
            plotSingularValues(arr_SingularValues, k_limits, t_limits);
            plotMinimumSingularValues(vMinimumSingularValues, k_limits, t_limits)
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
            plotResiduals(vMinimumResidual_QR, k_limits, t_limits);
            %plotResiduals(vMinimumResidual_SVD, k_limits);
        end
        
        metric = vMinimumResidual_QR;
        
end


LineBreakLarge();

% If only one subresultant exists, use an alternative method.
if (upperLimit_k - lowerLimit_k == 0 )
    
    % Get the metric to compute the degree of the GCD using only one
    % subresultant matrix.
    metric_OneSubresultants = arr_SingularValues{1};
    
    % Compute the degree of the GCD.
    t = Get_GCD_Degree_OneSubresultant_2Polys(metric_OneSubresultants);
    
    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    
    
    
else
    
    [t] = Get_GCD_Degree_MultipleSubresultants_2Polys(metric, k_limits);
    
    
    % Output all subresultants, all optimal alphas, all optimal thetas and all
    % geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
    alpha                   = vAlpha(t);
    theta                   = vTheta(t);
    GM_fx                   = vGM_fx(t);
    GM_gx                   = vGM_gx(t);
    
end




end






