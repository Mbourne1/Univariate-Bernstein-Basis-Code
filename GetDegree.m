function [t,out_subresultants_unprocessed,out_subresultants_preprocessed,out_alphas,out_thetas,out_gm_fxs,out_gm_gxs] = ...
    GetDegree(fx,gx)
% Get degree t of the AGCD d(x) of input polynomials f(x) and g(x)
%
%
%                       Inputs.
%
% fx : coefficients of polynomial f, expressed in Bernstein Basis
%
% gx : coefficients of polynomail g, expressed in Bernstein Basis
%
%                       Outputs.
%
% degree_calc - The calculated degree by various methods
%
% out_subresultants_unprocessed - All unprocesed subresultants S_{k} for k
% = 1,...,min(m,n)
%
% out_subresultants_preprocessed - All processed subresultants S_{k} for k
% = 1,...,min(m,n)
%
% out_alphas - All calculated optimal values of alpha from S_{k} for k =
% 1,...,min(m,n)
%
% out_thetas - All calculated optimal values of theta from S_{k} for k =
% 1,...,min(m,n)
%
% out_gm_fxs - All calculated geometric means of each C_{k}(f)
%
% out_gm_gxs - All calculated geometric means of each C_{k}(g)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global BOOL_PREPROC
global PLOT_GRAPHS
global THRESHOLD
global MIN_DELTA_MAG_ROW_SUMS

if isempty(BOOL_PREPROC) || isempty(PLOT_GRAPHS) ...
        isempty(THRESHOLD) || isempty(MIN_DELTA_MAG_ROW_SUMS)
    error('err')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get degree m of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get degree n of polynomial g(x)
[r,~] = size(gx);
n = r - 1;

fprintf('Begin : Get Degree m = % i , n = % i \n\n' ,m,n)


% get minimum degree of f(x) and g(x)
min_mn = min(m,n);

%% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
alpha_vec    =   zeros(min_mn,1);
theta_vec    =   zeros(min_mn,1);
gm_fx_vec    =   zeros(min_mn,1);
gm_gx_vec    =   zeros(min_mn,1);

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
vMinimumResidual             =   zeros(min_mn,1);

% Initialise a vector to store max/min diagonal entry in the upper
% triangular matrix R1_{k} from the QR decomposition of S_{k} for
% k=1,...,min(m,n)
vRatio_MaxMin_Diagonals_R    =   zeros(min_mn,1);

vMinimumSingularValues = zeros(min_mn,1);

% Initialise a vector to store max row sum / min row sum of the rows of
% R1_{k} from the QR decomposition of S_{k} for k=1,...,min(m,n)
vRatio_MaxMin_RowNorm_R  =   zeros(min_mn,1);

% Stores Data from QR decomposition.
Data_RowNorm    = [];
Data_DiagNorm   = [];

% Initialise arrays to store each subresultant matrix, both processed and
% unprocessed.
% Use Cell since the subresultants stored in the array are of varying
% sizes.
Sylvester_array_unproc          = cell(min_mn,1);
Sylvester_array_preprocessed    = cell(min_mn,1);

%% Loop

% For each subresultant $S_{k}$
for k = 1:1:min_mn
    
    % Get Unprocessed partitions (Including Geometric Mean)
    C_f_unproc = BuildT1(fx,n-k);
    C_g_unproc = BuildT1(gx,m-k);
    
    switch BOOL_PREPROC
        case 'y' % Include Preprocessors
            
            % Reason for performing BuildToeplitz before taking geometric
            % mean is that the MATLAB function geomean requires the
            % unprocessed matrix to calculate the GM.
            
            [lambda, mu] = GetGeometricMean(fx,gx,k);
            gm_fx_vec(k) = lambda;
            gm_gx_vec(k) = mu;
            
            % Divide normalised polynomials in Bernstein basis by geometric
            % means.
            fx_n    = fx / gm_fx_vec(k);
            gx_n    = gx / gm_gx_vec(k);
            
            % Divide the Sylvester Matrix partitions by Geometric mean.
            C_f_unproc_gm = C_f_unproc ./ gm_fx_vec(k);
            C_g_unproc_gm = C_g_unproc ./ gm_gx_vec(k);
            
            % Build subresultant S_{k}, and add to array of Sk
            Sylvester_array_unproc{k} = [C_f_unproc_gm C_g_unproc_gm];
            
            % For each coefficient ai of F, obtain the max and min such that F_max =
            % [max a0, max a1,...] and similarly for F_min, G_max, G_min
            [F_max,F_min,G_max,G_min] = GetMaxMin(Sylvester_array_unproc{k},m,n,k);
            
            % Calculate the optimal value of alpha and theta for the kth
            % subresultant matrix.
            [alpha_vec(k),theta_vec(k)] = OptimalAlphaTheta(F_max,F_min,G_max,G_min);
            
            % If linprog fails to find alpha and theta, and we have
            % assigned value 1. Best approximation is previous value of
            % alpha and theta when it exists.
            
        case 'n'
            % Exclude preprocessors
            
            % Dont normalise by geometric mean.
            fx_n = fx;
            gx_n = gx;
            
            % Dont obtain optimal alpha and theta.
            alpha_vec(k) = 1;
            theta_vec(k) = 1;
            gm_fx_vec(k) = 1;
            gm_gx_vec(k) = 1;
        otherwise
            error('bool_preproc must be either (y) or (n)')
    end
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = fx_n.*(theta_vec(k).^(0:1:m)');
    gw = gx_n.*(theta_vec(k).^(0:1:n)');
    
    Sk  =   BuildSubresultant(fw,gw,k,alpha_vec(k));
    
    % Add Sk to the array of preprocessed Sk
    Sylvester_array_preprocessed{k} = Sk;
    
    % Using
    vMinimumSingularValues(k) = min(svd(Sk));
    
    % Get QR Decomposition of the sylvester matrix S_{k}(f(\theta,w),g(\theta,w))
    [~,R] = qr(Sk);
    
    % Take absolute values of R_{k}
    R = abs(R);
    
    % Get number of rows in R1_{k}
    [nRowsR1,~] = size(diag(R));
    
    % Obtain R1 the top square of the R matrix.
    R1 = R(1:nRowsR1,1:nRowsR1);
    
    % Get Norms of each row in the matrix R1
    vR1_RowNorms = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    vR1_DiagNorm = diag(R1)./norm(diag(R1));
    
    % Scatter Plot Data
    ks = k.*ones(size(vR1_RowNorms));
    ns = 1:1:size(vR1_RowNorms,1);
    
    % Form a triple of [ks, the value of QR_RowNorm, and the index of the value of
    % the row of R1 corresponding to QR_RowNorm].
    
    X = [ks vR1_RowNorms ns'];
    Data_RowNorm = [Data_RowNorm; X];
    
    X2 = [ks vR1_DiagNorm ns'];
    Data_DiagNorm = [Data_DiagNorm;X2];
    
    % Get max:min diag elem of S_{k}
    vRatio_MaxMin_Diagonals_R(k) = max(diag(R1))./min(diag(R1));
    vRatio_MaxMin_Diagonals_R(k) = sanitize(vRatio_MaxMin_Diagonals_R(k));
    
    % Get max:min rownorm r_{i}/r_{j} of S_{k}
    vRatio_MaxMin_RowNorm_R(k) = max(vR1_RowNorms)./min(vR1_RowNorms);
    vRatio_MaxMin_RowNorm_R(k) = sanitize(vRatio_MaxMin_RowNorm_R(k));
    
    % Get the minimal residuals for each subresultant S_{k}.
    vMinimumResidual(k) = GetMinimalDistance(Sk);
    
end

%% Analyse Max:Min Row Norms for each subresultant

% Get the change in the ratios from one subresultant to the next.
vDelta_MaxMin_RowNorm_R = abs(diff(log10(vRatio_MaxMin_RowNorm_R)));

% Get the maximum change in rowsum ratio and its index
[max_Delta_MaxMin_RowSum_R,indexMaxChange_RowNorm] = max(vDelta_MaxMin_RowNorm_R);

%% Anaysis of Max:Min Diagonal entries of each subresultant

% Get the change in the ratios of diagonal elements from one subresultant
% to the next.
vDelta_MaxMin_Diag_R = abs(diff(log10(vRatio_MaxMin_Diagonals_R)));

% Get the maximum change in diag ratio and its index
[max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = max(vDelta_MaxMin_Diag_R);

if min(m,n) == 1
    t = Get_Rank_One_Subresultant(Sk);
    
    out_subresultants_unprocessed = Sylvester_array_unproc{1};
    out_subresultants_preprocessed = Sylvester_array_preprocessed{1};
    out_alphas = alpha_vec;
    out_thetas = theta_vec;
    out_gm_fxs = gm_fx_vec;
    out_gm_gxs = gm_gx_vec;
    return;
end


% Set a condition for which we consider the maximum change in row sums to
% significant or insignificant.
if abs(max_Delta_MaxMin_Diag_R) < THRESHOLD
    
    % Check to see if all subresultants are rank deficient in which case
    % the degree of the GCD is min(m,n)
    
    if mean(log10(vRatio_MaxMin_Diagonals_R)) < THRESHOLD
        
        %all subresultants are full rank
        t = 0;
        
        out_subresultants_unprocessed = 0;
        out_subresultants_preprocessed = 0;
        out_alphas = 0;
        out_thetas = 0;
        out_gm_fxs = 0;
        out_gm_gxs = 0;
        fprintf('All subresultants are Non-Singular \n')
        return;
    else
        % All subresultants are rank deficient
        fprintf('All subresultants are Singular \n')
        t = min(m,n);
    end

else % The degree of the GCD is somewhere between 1 and min(m,n)
    
    % Get degree of GCD by ration of max:min row norms
    degree_calc_1 = indexMaxChange_RowNorm;
    
    % Get degree by ratio of max/min element in N, where N is a vector of the
    % norms of each row of R1.
    degree_calc_2 = indexMaxChange_Ratio_DiagsR;
    
    % Get degree by minimal residual obtained by removing each column from each
    % subresultant, giving Ak x = ck. When residual is 'close to zero'. Matrix is rank
    % deficient. (Note if zero, log10(0) = inf. so use fudge factor)
    A = log10(vMinimumResidual);
    [~,degree_calc_3] = max(abs(diff(A)));
    
    % Take the mode of the three degree calculation methods
    t = mode([degree_calc_1,degree_calc_2,degree_calc_3]);
    
    
    
    
end

Plotting();
%%


%%
% Outputs

% Output just corresponding to calculated value of the degree.
% Output subresultant S_{t}, alpha_{t}, theta_{t}, and corresponding
% geometric means.
% out_subresultant_unprocessed = cell2mat(Sylvester_array_unproc(deg_calc));
% out_subresultant_preprocessed = cell2mat(Sylvester_array_preprocessed(deg_calc));
% out_alpha = alpha_vec(deg_calc);
% out_theta = theta_vec(deg_calc);
% out_gm_fx = gm_fx_vec(deg_calc);
% out_gm_gx = gm_gx_vec(deg_calc);

% Output all subresultants, all optimal alphas, all optimal thetas and all
% geometric means for each subresultant S_{k} where k = 1,...,min(m,n)
out_subresultants_unprocessed = Sylvester_array_unproc;
out_subresultants_preprocessed = Sylvester_array_preprocessed;
out_alphas = alpha_vec;
out_thetas = theta_vec;
out_gm_fxs = gm_fx_vec;
out_gm_gxs = gm_gx_vec;


end