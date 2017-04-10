function [t, alpha, beta, theta, GM_fx, GM_gx, GM_hx] = ...
    Get_GCD_Degree_3Polys(fx, gx, hx, limits_t)
% GetGCD_Degree_2Polys(fx,gx)
%
% Get degree t of the AGCD d(x) of input polynomials f(x) and g(x)
%
%
% % Inputs.
%
% fx : (Vector) coefficients of polynomial f(x)
%
% gx : (Vector) coefficients of polynomail g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% limits_t : [Int Int] Set the upper and lower bound of the degree of the
% GCD of polynomials f(x) and g(x). Usually used when using o_roots() where
% prior information is known about the upper and lower bound due to number
% of distinct roots.
%
%
% % Outputs.
%
% t : (Int) Degree of GCD of f(x) and g(x)
%
% alpha : (Float) optimal value of alpha
%
% beta : (Float) 
%
% theta : (Float) optimal value of theta
%
% GM_fx : (Float) Geometric mean of entries of f(x)
% 
% GM_gx : (Float) Geometric mean of entries of g(x)
%
% GM_hx : (Float) Geometric mean of entries of h(x)
%

addpath 'Preprocessing'
addpath 'Sylvester Matrix'

% Get degree of polynomail f(x)
m = GetDegree(fx);
n = GetDegree(gx);
o = GetDegree(hx);

% Set upper and lower limit for computing the degree of the GCD. Note may
% be best to set to degree limits or may be best to set to 1 & min(m,n)
lowerLimit_k = 1;
upperLimit_k = min([m,n,o]);
limits_k = [lowerLimit_k upperLimit_k];

% Get the number of subresultants which must be constructed.
nSubresultants = upperLimit_k - lowerLimit_k +1 ;

% %
% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha = zeros(nSubresultants, 1);
vBeta = zeros(nSubresultants, 1);
vTheta = zeros(nSubresultants, 1);
vGM_fx = zeros(nSubresultants, 1);
vGM_gx = zeros(nSubresultants, 1);
vGM_hx = zeros(nSubresultants, 1);

arr_Sk = cell(nSubresultants, 1);
arr_R1 = cell(nSubresultants, 1);

vMaxDiagR1 = zeros(nSubresultants,1);
vMinDiagR1 = zeros(nSubresultants,1);
vMaxRowNormR1 = zeros(nSubresultants,1);
vMinRowNormR1 = zeros(nSubresultants,1);

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
vMinimumResidual_QR  = zeros(nSubresultants,1);
vMinimumResidual_SVD = zeros(nSubresultants,1);

% Initialise a vector to store the minimum singular values \sigma_{i} of each
% subresultant S_{i}(f,g), where i is between the upper and lower bound.
vMinimumSingularValues = zeros(nSubresultants,1);




% For each subresultant $S_{k}$
for k = lowerLimit_k:1:upperLimit_k
    
    i = k - lowerLimit_k + 1;
    
    [vGM_fx(i), vGM_gx(i), vGM_hx(i), alpha, beta, theta ] = Preprocess_3Polys(fx, gx, hx, k);
    
    vAlpha(i) = alpha;
    vBeta(i) = beta; 
    vTheta(i) = theta;
    
    display(alpha)
    display(beta) 
    display(theta)
    
    % 18/04/2016
    % Given the previous geometric mean of f(x) calculate the new geometric
    % mean by my new method
    if i>1
        %GM_fx = GetGeometricMeanFromPrevious(fx , vGM_fx(i-1) , m , n-k);
        %GM_gx = GetGeometricMeanFromPrevious(gx , vGM_gx(i-1) , n , m-k);
    end
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx./ vGM_fx(i);
    gx_n = gx./ vGM_gx(i);
    hx_n = hx./ vGM_hx(i);
    
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = GetWithThetas(fx_n, vTheta(i));
    gw = GetWithThetas(gx_n, vTheta(i));
    hw = GetWithThetas(hx_n, vTheta(i));
    
    % Build the k-th subresultant   
    % Note : We dont include any alpha yet 
    arr_Sk{i} = BuildSubresultant_3Polys(fw, alpha.*gw, beta.*hw,k);
    
    
    
    % Get the matrix R1 from the QR Decomposition of S
    arr_R1{i} = GetR1(arr_Sk{i});
    
      
    
    
end % End of for

% %
% %
% %
% Choose a metric to determine the degree of the GCD.
% R1 Row Norms
% R1 Row Diagonals
% Singular Values
% Residuals

global SETTINGS
fprintf(sprintf('Metric used to compute degree of GCD : %s \n', SETTINGS.RANK_REVEALING_METRIC));

switch SETTINGS.RANK_REVEALING_METRIC
    case 'Row Norms'
        
        for i = 1:1:nSubresultants
            
            % Get Norms of each row in the matrix R1
            vR1_RowNorms = sqrt(sum(arr_R1{i}.^2,2))./norm(R1);

            % Get maximum and minimum row norms of rows of R1.
            vMaxRowNormR1(i) = max(vR1_RowNorms);
            vMinRowNormR1(i) = min(vR1_RowNorms);
            
        end
        
        metric = vMaxRowNormR1./vMinRowNormR1;
        
        plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t);

    case 'Row Diagonals'
        
        for i = 1:1:nSubresultants
            
            % Get Norms of diagonal eleements of R1
            vDiagsR1 = diag(arr_R1{i});

            % Get maximum and minimum row diagonals of R1
            vMaxDiagR1(i) = max(vDiagsR1);
            vMinDiagR1(i) = min(vDiags1);
            
        end
        
        plotMaxMinRowDiagonals(vMaxDiagR1, vMinDiagR1, limits_k, limits_t);
        
        metric = vMaxDiagR1./vMinDiagR1;
        
    case 'Singular Values'
       
        arr_SingularValues = cell(nSubresultants,1);
        
        for i = 1 : 1 : nSubresultants
        
            % Get singular values of S_{k}
            arr_SingularValues{i} = svd(arr_Sk{i});

            % Get the minimal singular value from S_{k}
            vMinimumSingularValues(i) = min(arr_SingularValues{i});
        
        end
        
        metric = vMinimumSingularValues;
        
        plotSingularValues(arr_SingularValues, limits_k, limits_t);
        plotMinimumSingularValues(vMinimumSingularValues, limits_k, limits_t);
        
    case 'Residuals'
        % Get the minimal residuals for each subresultant S_{k}.
        vMinimumResidual_QR(i) = GetMinimalDistance(arr_Sk{i},'QR');
        vMinimumResidual_SVD(i) = GetMinimalDistance(arr_Sk{i},'SVD');
        
        error('Code not developed')
        
    otherwise
        error('err')
end


% % Analysis of Minimum Singular values

if (upperLimit_k == lowerLimit_k)
    
    % Set \alpha and \theta
    alpha = vAlpha(1);
    beta = vBeta(1);
    theta = vTheta(1);
    
    % Set Geometric means
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    GM_hx = vGM_hx(1);
    
    % Set degree of GCD
    t = upperLimit_k;
    
    return;
end


% If only one subresultant exists, use an alternative method.
if (upperLimit_k - lowerLimit_k == 0 )
    
        % Set the degree of the GCD
        t = Get_GCD_Degree_OneSubresultant(vSingularValues);
    
        % Set \alpha and \theta
        alpha = vAlpha(1);
        beta = vBeta(1);
        theta = vTheta(1);
        
        % Set geometric means
        GM_fx = vGM_fx(1);
        GM_gx = vGM_gx(1);
        GM_hx = vGM_hx(1);
        
        return;
        
    
    
else
    
    [t] = Get_GCD_Degree_MultipleSubresultants_2Polys(metric, limits_k);
    
    
    % % Graph Plotting
    %PlotGraphs_GCDDegree()
    
    % %
    % Outputs
    
    % Output all subresultants, all optimal alphas, all optimal thetas and all
    % geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
    
    alpha = vAlpha(t);
    theta = vTheta(t);
    
    GM_fx = vGM_fx(t);
    GM_gx = vGM_gx(t);
    GM_hx = vGM_hx(t);
    
end




end


function data = AddToResults(data,vector,k)
% Given the vector of data, it is necessary to include in a 3 columned
% matrix, along with the corresponding value k to indicate that the value
% came from the k-th subresultant, and an index 1:1:nEntries.


% Get Number of entries in vector
nEntries = size(vector,1);

% Get a vector of k
v_k = k.* ones(nEntries,1);

% Get a vector of 1:1:n
vec_n = 1:1:nEntries;

% Form a triple of [ks, the value of QR_RowNorm, and the index of the value of
% the row of R1 corresponding to QR_RowNorm].

temp_data = [v_k vector vec_n'];
data = [data; temp_data];

end




