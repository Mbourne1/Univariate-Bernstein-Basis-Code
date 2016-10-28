function [t,alpha, theta,GM_fx,GM_gx] = ...
    GetGCD_DegreeByNewMethod(fx,gx)
% Get degree of the AGCD of input polynomials f(x) and g(x)
%
% Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% Outputs.
%
% deg_calc - The calculated degree by various methods
%
% out_subresultants_unprocessed - All unprocesed subresultants S_{k} for k
% = 1,...,min(m,n)
%
% out_subresultants_preprocessed - All processed subresultants S_{k} for k
% = 1,...,min(m,n)
%
% out_alphas : All calculated optimal values of alpha from S_{k} for k =
% 1,...,min(m,n)
%
% out_thetas : All calculated optimal values of theta from S_{k} for k =
% 1,...,min(m,n)
%
% out_gm_fxs : All calculated geometric means of each C_{k}(f)
%
% out_gm_gxs : All calculated geometric means of each C_{k}(g)
%
%



% Get degree m of polynomial f
m = GetDegree(fx);

% Get degree n of polynomial g
n = GetDegree(gx);

% get minimum degree of f and g
min_mn = min(m,n);

%
lower_lim = 1;
upper_lim = min_mn;


% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(1,min_mn);
vTheta    =   zeros(1,min_mn);
vGM_fx    =   zeros(1,min_mn);
vGM_gx    =   zeros(1,min_mn);


% Initialise vectors to store values calculated from each subresultant
% S_{k} for k = 1,...,min(m,n).

vMaxDiagR1  = zeros(1,min_mn);
vMinDiagR1  = zeros(1,min_mn);
vMaxRowNormR1   = zeros(1,min_mn);
vMinRowNormR1   = zeros(1,min_mn);
vMinimumResidual    = zeros(1,min_mn);
vMinimumSingularValues  = zeros(1,min_mn);

% Stores Data from QR decomposition.
Data_RowNorm = [];
Data_DiagNorm = [];


% For each subresultant $$S_{k}$$
for k = 1:1:min_mn
    
    % if this is the first subresultant, build from scratch
    if (k == 1)
        % Build the partitions of S_{k}(f,g)
        DT1Q1_unproc = BuildDT1Q1(fx,n-k);
        DT2Q2_unproc = BuildDT1Q1(gx,m-k);
    else
        % Get unprocessed partitions
        DT1Q1_unproc = BuildDT1Q1_fromPrev(m,n,k-1,DT1Q1_unproc);
        DT2Q2_unproc = BuildDT1Q1_fromPrev(n,m,k-1,DT2Q2_unproc);
        
    end
    
    % Preprocess polynomials f(x) and g(x) in the Sylvester matrix
    % S_{k}(f,g)
    [vGM_fx(k), vGM_gx(k), vAlpha(k), vTheta(k)] = Preprocess(fx,gx,k);
    
    % Divide the Sylvester Matrix partitions by Geometric mean.
    DT1Q1_unproc = DT1Q1_unproc ./ vGM_fx(k);
    DT2Q2_unproc = DT2Q2_unproc ./ vGM_gx(k);
    
    % Divide normalised polynomials in Bernstein basis by
    % geometric means.
    fx_n = fx./vGM_fx(k);
    gx_n = gx./vGM_gx(k);
    
    % Calculate the coefficients of the modified Bernstein basis
    % polynomials F2 and G2. Multiply G2 by alpha.
    fw_n = GetWithThetas(fx_n,vTheta(k));
    gw_n = GetWithThetas(gx_n,vTheta(k));
    
    % Construct the kth subresultant matrix for the optimal values of alpha
    % and theta.
    if (k == 1)
        % first subresultant must be constructed in the conventional sense.
        DT1Q1 = BuildDT1Q1(fw_n,n-k);
        DT2Q2 = BuildDT1Q1(gw_n,m-k);
        
        % Build subresultant
        DTQ = [DT1Q1 vAlpha(k).*DT2Q2]   ;
        
    else
        % Subsequent subresultants can be built using my build up method.
        DT1Q1 = BuildDT1Q1_fromPrev(m,n,k-1,DT1Q1);
        DT2Q2 = BuildDT1Q1_fromPrev(n,m,k-1,DT2Q2);
        DTQ = [DT1Q1 DT2Q2];
        
    end
    
    
    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(DTQ);
    
    % Take absolute values.
    R_abs = abs(R);
    
    % Get number of rows in R1
    [nRowsR1,~] = size(diag(R_abs));
    
    % Obtain R1, the top square of the R matrix
    R1_abs = R_abs(1:nRowsR1,1:nRowsR1);
    
    % Get Norms of each row in matrix R1.
    R1_RowNorm = sqrt(sum(R1_abs.^2,2))./norm(R1_abs);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1_abs)./norm(diag(R1_abs));
    
    %Scatter Plot Data
    ks = k.*ones(size(R1_RowNorm));
    ns = 1:1:size(R1_RowNorm,1);
    
    % Form a triple of [ks, the value of QR_RowNorm, and the index of the
    % value of the row of R1 corresponding to QR_RowNorm]. EG.
    %  [1   0.015  1
    %   1   0.156  2 2 ...]
    X = [ks R1_RowNorm ns'];
    Data_RowNorm = [Data_RowNorm; X];
    
    X2 = [ks R1_DiagNorm ns'];
    Data_DiagNorm = [Data_DiagNorm;X2];
    
    %
    % Get ratio of max diag elem of S_{k} to min diag elemente of S_{k}
    
    vMaxDiagR1(k) = max(diag(R1_abs));
    vMinDiagR1(k) = min(diag(R1_abs));
    
    vMaxRowNormR1(k) = max(R1_RowNorm);
    vMinRowNormR1(k) = min(R1_RowNorm);
    
    
    vMinimumResidual(k) = GetMinimalDistance(DTQ);
    
    vSingularValues = svd(DTQ);
    
    vMinimumSingularValues(k) = min(vSingularValues);
    
end % End of for

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



% Check to see if only one subresultant exists, ie if m or n is equal
% to one

if min_mn == 1
    t = GetGCDDegree_OneSubresultant(vSingularValues);
    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    return;
    
else
    
    %
    %
    [t] = GetGCDDegree_MultipleSubresultants(metric,lower_lim);
       
    % Outputs
    
    % Output just corresponding to calculated value of the degree. Output
    % subresultant S_{t}, alpha_{t}, theta_{t}, and corresponding geometric
    % means.
    alpha = vAlpha(t);
    theta = vTheta(t);
    GM_fx = vGM_fx(t);
    GM_gx = vGM_gx(t);
    
    % % Graph Plotting
    PlotGraphs()
    
end








end


