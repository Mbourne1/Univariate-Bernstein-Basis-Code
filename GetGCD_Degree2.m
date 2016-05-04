function [t,alpha, theta,GM_fx,GM_gx] = ...
    GetGCD_Degree2(fx,gx,deg_limits)
% GetGCD_Degree(fx,gx)
%
% Get degree t of the AGCD d(x) of input polynomials f(x) and g(x)
%
%
% % Inputs.
%
% fx : coefficients of polynomial f, expressed in Bernstein Basis
%
% gx : coefficients of polynomail g, expressed in Bernstein Basis
%
% deg_limits : set the upper and lower bound of the degree of the
% GCD of polynomials f(x) and g(x). Usually used when using o_roots() where
% prior information is known about the upper and lower bound due to number
% of distinct roots. 
%
%
% % Outputs.
%
% t : Degree of GCD of f(x) and g(x)
%
% alpha : optimal value of alpha
%
% theta : optimal value of theta
%

% Get global variables

global PLOT_GRAPHS
global THRESHOLD


% Get degree m of polynomial f(x)
m = GetDegree(fx);

% Get degree n of polynomial g(x)
n = GetDegree(gx);

% If the number of distinct roots in f(x) is one, then the degree of the
% GCD of f(x) and f'(x) = m-1 = n.
lower_lim = deg_limits(1);
upper_lim = deg_limits(2);

if (upper_lim == lower_lim)
    t = upper_lim;
    return;
end

% if the lower limit is zero, ie - deg(GCD) = 0, then set to 1, since we
% are only interested in Sylvester matrices S_{1},...
if lower_lim == 0
    lower_lim = 0;
else 
end

% Get the number of subresultants which must be constructed.
nSubresultants = upper_lim - lower_lim +1 ;

% %
% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(nSubresultants,1);
vTheta    =   zeros(nSubresultants,1);
vGM_fx    =   zeros(nSubresultants,1);
vGM_gx    =   zeros(nSubresultants,1);

vMax_Diag_R1 = zeros(nSubresultants,1);
vMin_Diag_R1 = zeros(nSubresultants,1);
vMaxRowNormR1 = zeros(nSubresultants,1);
vMinRowNormR1 = zeros(nSubresultants,1);

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
vMinimumResidual             =   zeros(nSubresultants,1);

% Initialise a vector to store max/min diagonal entry in the upper
% triangular matrix R1_{k} from the QR decomposition of S_{k} for
% k=1,...,min(m,n)
vRatio_MaxMin_Diagonals_R    =   zeros(nSubresultants,1);


vMinimumSingularValues = zeros(nSubresultants,1);

% Stores Data from QR decomposition.
Data_RowNorm    = [];
Data_DiagNorm   = [];


%% Loop

% For each subresultant $S_{k}$
for k = lower_lim:1:upper_lim
    
    i = k - lower_lim + 1;
    
    [vGM_fx(i), vGM_gx(i),vAlpha(i),vTheta(i)] = Preprocess(fx,gx,k);
    
    % 18/04/2016
    % Given the previous geometric mean of f(x) calculate the new geometric
    % mean by my new method
    if i>1
        GM = GetGeometricMeanFromPrevious(fx , vGM_fx(i-1) , m , n-k);
    end
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx./vGM_fx(i);
    gx_n = gx./vGM_gx(i);
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = GetWithThetas(fx_n,vTheta(i));
    gw = GetWithThetas(gx_n,vTheta(i));
    
    % Build the k-th subresultant
    Sk = BuildSubresultant(fw,vAlpha(i).*gw,k);
    
    % Get the minimal singular value from S_{k}
    vMinimumSingularValues(i) = min(svd(Sk));
    
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
    vMax_Diag_R1(i) = max(vDiagsR1);
    vMin_Diag_R1(i) = min(vDiagsR1);
    
    % Get maximum and minimum row norms of rows of R1.
    vMaxRowNormR1(i) = max(vR1_RowNorms);
    vMinRowNormR1(i) = min(vR1_RowNorms);
    
    % Get the minimal residuals for each subresultant S_{k}.
    vMinimumResidual(i) = GetMinimalDistance(Sk);
    
    
end

% End of Loop

% Get ratio of max:min Row Norms of R1
vRatio_MaxMin_RowNorm_R = vMaxRowNormR1./ vMinRowNormR1;
vRatio_MaxMin_RowNorm_R = sanitize(vRatio_MaxMin_RowNorm_R);

% Get ratio of max:min diagonals of R1
vRatio_MaxMin_Diagonals_R = vMax_Diag_R1 ./ vMin_Diag_R1;
vRatio_MaxMin_Diagonals_R = sanitize(vRatio_MaxMin_Diagonals_R);

% % Analyse Max:Min Row Norms for each subresultant
% Get the change in the ratios from one subresultant to the next.
vDelta_MaxMin_RowNorm_R = abs(diff(log10(vRatio_MaxMin_RowNorm_R)));

% Get the maximum change in rowsum ratio and its index
[max_Delta_MaxMin_RowSum_R,indexMaxChange_RowNorm] = max(vDelta_MaxMin_RowNorm_R);

% % Anaysis of Max:Min Diagonal entries of each subresultant

% Get the change in the ratios of diagonal elements from one subresultant
% to the next.
vDelta_MaxMin_Diag_R = abs(diff(log10(vRatio_MaxMin_Diagonals_R)));

% Get the maximum change in diag ratio and its index
[max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = max(vDelta_MaxMin_Diag_R);

% % Analysis of Minimum Singular values


% If only one subresultant exists, use an alternative method.
if (upper_lim -lower_lim == 0 ) % If only one Subresultant Exists
     t = GetRank_One_Subresultant(Sk);
     alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    return;
end

[ProblemType,t] = GetProblemType(vMinimumSingularValues,lower_lim);


%% Graph Plotting
plotgraphs()


% %  
% Outputs

% Output all subresultants, all optimal alphas, all optimal thetas and all
% geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
alpha                   = vAlpha(t-lower_lim+1);
theta                   = vTheta(t-lower_lim+1);
GM_fx                   = vGM_fx(t-lower_lim+1);
GM_gx                   = vGM_gx(t-lower_lim+1);
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

function R1 = GetR1(Sk)
% Get the square upper triangular matrix R1 from the QR Decomposition of
% the kth subresultant matrix Sk
%
% Inputs.
%
% Sk : k-th Subresultant matrix S_{k}(f,g)

% Get QR Decomposition of the sylvester matrix S_{k}(f(\theta,w),g(\theta,w))
[~,R] = qr(Sk);

% Take absolute values of R_{k}
R = abs(R);

% Get number of rows in R1_{k}
[nRowsR1,~] = size(diag(R));

% Obtain R1 the top square of the R matrix.
R1 = R(1:nRowsR1,1:nRowsR1);


end


function [ProblemType , t] = GetProblemType(vMinimumSingularValues,lower_lim)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.



global THRESHOLD

min_mn = lower_lim + length(vMinimumSingularValues) - 1;

% Get the differences between minimum singular value S_{k} and S_{k+1}
vDeltaMinSingVal = abs(diff(log10(vMinimumSingularValues)));
vDeltaMinSingVal(isnan(vDeltaMinSingVal)) = 0 ;
vDeltaMinSingVal(vDeltaMinSingVal==inf)=0;

% Get the index of the largest change in minimum singular values
[maxChangeSingularValues, indexMaxChange] = max(vDeltaMinSingVal);

if  abs(maxChangeSingularValues) < THRESHOLD
    
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = log10(norm(vMinimumSingularValues));
    
    if  avgMinSingularValue < THRESHOLD
        % If all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        ProblemType = 'NonSingular';
        fprintf('All Subresultants are NonSingular \n')
        t = min_mn;
    else
        % if all singular values are not close to zero, then full rank, degree
        % of gcd is 0
        ProblemType = 'Singular';
        fprintf('All subresultants are Singular \n')
        t = 0;
    end
else
    % maxChange is signifcant
    ProblemType = 'Mixed';
    fprintf('Mixed \n')
    t = lower_lim + indexMaxChange - 1;
end

end
