function [t,alpha, theta,GM_fx,GM_gx] = ...
    GetGCD_Degree(fx,gx)
% GetGCD_Degree(fx,gx)
%
% Get degree t of the AGCD d(x) of input polynomials f(x) and g(x)
%
%
% Inputs.
%
% fx : coefficients of polynomial f, expressed in Bernstein Basis
%
% gx : coefficients of polynomail g, expressed in Bernstein Basis
%
%                       Outputs.
%
% t : Degree of GCD of f(x) and g(x)
%
% alpha : optimal value of alpha
%
% theta : optimal value of theta
%

% Get global variables

% Get degree m of polynomial f(x)
m = GetDegree(fx);

% Get degree n of polynomial g(x)
n = GetDegree(gx);

% get minimum degree of f(x) and g(x)
min_mn = min(m,n);

% Set upper and lower limit
lower_lim = 1;
upper_lim = min_mn;


%% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(min_mn,1);
vTheta    =   zeros(min_mn,1);
vGM_fx    =   zeros(min_mn,1);
vGM_gx    =   zeros(min_mn,1);

vMaxDiagR1 = zeros(min_mn,1);
vMinDiagR1 = zeros(min_mn,1);
vMaxRowNormR1 = zeros(min_mn,1);
vMinRowNormR1 = zeros(min_mn,1);

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
vMinimumResidual             =   zeros(min_mn,1);

% Initialise a vector to store minimal singular values of each S_{k}
vMinimumSingularValues = zeros(min_mn,1);

% Stores Data from QR decomposition.
Data_RowNorm    = [];
Data_DiagNorm   = [];

%%

% For each subresultant $S_{k}$
for k = 1:1:min_mn
    
    
    [vGM_fx(k), vGM_gx(k),vAlpha(k),vTheta(k)] = Preprocess(fx,gx,k);
    
    % 18/04/2016
    % Given the previous geometric mean of f(x) calculate the new geometric
    % mean by my new method
    if k>1
        GM_fx = GetGeometricMeanFromPrevious(fx , vGM_fx(k-1) , m , n-k);
        GM_gx = GetGeometricMeanFromPrevious(gx , vGM_gx(k-1) , n , m-k);
    end
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx./vGM_fx(k);
    gx_n = gx./vGM_gx(k);
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = GetWithThetas(fx_n,vTheta(k));
    gw = GetWithThetas(gx_n,vTheta(k));
    
    % Build the k-th subresultant
    Sk = BuildSubresultant(fw,vAlpha(k).*gw,k);
    
    % Get the singular values of S_{k}
    vSingularValues = svd(Sk);
    
    % Get the minimal singular value from S_{k}
    vMinimumSingularValues(k) = min(vSingularValues);
    
    % Get the matrix R1 from the QR Decomposition of S
    R1 = GetR1(Sk);
    
    % Get Norms of each row in the matrix R1
    vR1_RowNorms = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get Norms of diagonal eleements of R1
    vDiagsR1 = diag(R1);
    vDiagsR1_norm = vDiagsR1 ./ norm(vDiagsR1);
    
    % Insert Row Norm data into a matrix
    Data_RowNorm = AddToResults(Data_RowNorm,vR1_RowNorms,k);
    
    % Insert diagonals data into a matrix
    Data_DiagNorm = AddToResults(Data_DiagNorm,vDiagsR1_norm,k);
    
    % Get maximum and minimum row diagonals of R1
    vMaxDiagR1(k) = max(vDiagsR1);
    vMinDiagR1(k) = min(vDiagsR1);
    
    % Get maximum and minimum row norms of rows of R1.
    vMaxRowNormR1(k) = max(vR1_RowNorms);
    vMinRowNormR1(k) = min(vR1_RowNorms);
    
    % Get the minimal residuals for each subresultant S_{k}.
    vMinimumResidual(k) = GetMinimalDistance(Sk);
    
    
end % End of For

% %
% %
% %
% Choose a metric to determine the degree of the GCD.
global SETTINGS
switch SETTINGS.METRIC
    case 'Row Norms'
        metric = vMaxRowNormR1./vMinRowNormR1;
        
    case 'Row Diagonals'
        metric = vMaxDiagR1./vMinDiagR1;
        
    case 'Singular Values'
        metric = vMinimumSingularValues;
end



%
% If only one subresultant exists, get GCD by alternative method
if min_mn == 1
    
    fprintf([mfilename ' : ' 'Only One Subresultant \n'])
    t = GetGCDDegree_OneSubresultant(vSingularValues);
    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    return;
    
else
    t = GetGCDDegree_MultipleSubresultants(metric,lower_lim);
    
    % % Graph Plotting
    PlotGraphs()
    % Outputs
    
    % Output all subresultants, all optimal alphas, all optimal thetas and all
    % geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
    
    alpha                   = vAlpha(t);
    theta                   = vTheta(t);
    GM_fx                   = vGM_fx(t);
    GM_gx                   = vGM_gx(t);
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




