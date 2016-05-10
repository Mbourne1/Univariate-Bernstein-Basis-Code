function [t,alpha, theta,GM_fx,GM_gx] = ...
    GetGCD_DegreeByNewMethod(fx,gx,deg_limits)
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


% Get lower limit of degree of gcd
lower_lim = deg_limits(1);
upper_lim = deg_limits(2);

% Get the number of subresultants to be built
nSubresultants = upper_lim - lower_lim +1 ;

% Get degree m of polynomial f
m = GetDegree(fx);

% Get degree n of polynomial g
n = GetDegree(gx);

% get minimum degree of f and g
min_mn = min(m,n);

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(1,nSubresultants);
vTheta    =   zeros(1,nSubresultants);
vGM_fx    =   zeros(1,nSubresultants);
vGM_gx    =   zeros(1,nSubresultants);


% Initialise vectors to store values calculated from each subresultant
% S_{k} for k = 1,...,min(m,n).
vMaxDiagR1= zeros(1,nSubresultants);
vMinDiagR1= zeros(1,nSubresultants);
vMaxRowNormR1= zeros(1,nSubresultants);
vMinRowNormR1 = zeros(1,nSubresultants);
vMinimumResidual = zeros(1,nSubresultants);
vMinimumSingularValues = zeros(1,nSubresultants);
% Stores Data from QR decomposition.
Data_RowNorm = [];
Data_DiagNorm = [];


% For each subresultant $$S_{k}$$
for k = lower_lim:1:upper_lim
    
    % Get index i, used for storing in vectors.
    i = k - lower_lim + 1;
    
    % If S_{k} is the first subresultant, then it must be built from
    % scratch
    if (i == 1)
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
    [vGM_fx(i), vGM_gx(i), vAlpha(i), vTheta(i)] = Preprocess(fx,gx,k);
    
    % Divide the Sylvester Matrix partitions by Geometric mean.
    DT1Q1_unproc = DT1Q1_unproc ./ vGM_fx(i);
    DT2Q2_unproc = DT2Q2_unproc ./ vGM_gx(i);
    
    % Divide normalised polynomials in Bernstein basis by
    % geometric means.
    fx_n = fx./vGM_fx(i);
    gx_n = gx./vGM_gx(i);
    
    % Calculate the coefficients of the modified Bernstein basis
    % polynomials F2 and G2. Multiply G2 by alpha.
    fw_n = GetWithThetas(fx_n,vTheta(i));
    gw_n = GetWithThetas(gx_n,vTheta(i));
    
    % Construct the kth subresultant matrix for the optimal values of alpha
    % and theta.
    if (i == 1)
        % first subresultant must be constructed in the conventional sense.
        DT1Q1 = BuildDT1Q1(fw_n,n-k);
        DT2Q2 = BuildDT1Q1(gw_n,m-k);
        
        % Build subresultant
        DTQ = [DT1Q1 vAlpha(i).*DT2Q2]   ;
        
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
    
    vMaxDiagR1(i) = max(diag(R1_abs));
    vMinDiagR1(i) = min(diag(R1_abs));
    
    vMaxRowNormR1(i) = max(R1_RowNorm);
    vMinRowNormR1(i) = min(R1_RowNorm);
    
    % For each subresultant Sk - Remove each column c_{k,i} in turn, obtain
    % the residuals c_{k}-A_{k}.
    
    % residualQR_vector :- stores all residuals for 1,...,m+n-k+2 for a
    % given subresultant
    
    % minResQR_vector :- stores only one residual (the min) for each
    % subresultant k=1,...min(m,n)
    
    vMinimumResidual(i) = GetMinimalDistance(DTQ);
    vMinimumSingularValues(i) = min(svd(DTQ));
    
end

% %
% %
% %
% %
[max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = Analysis(vMaxDiagR1 ./ vMinDiagR1);

% %
% %
% %
% %
[max_delta_mag_rowsum,indexMaxChange_RowNorm] = Analysis(vMaxRowNormR1 ./ vMinRowNormR1);


% %
% %
% %
% %
[max_delta_min_residuals,indexMaxChange_Residuals] = Analysis(vMinimumResidual);




% Check to see if only one subresultant exists, ie if m or n is equal
% to one

if nSubresultants == 1
    t = GetRankOneSubresultant(DTQ);
    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
    return;
end

%
%
[t] = GetProblemType(vMinimumSingularValues,lower_lim);

PlotGraphs();

% Outputs

% Output just corresponding to calculated value of the degree. Output
% subresultant S_{t}, alpha_{t}, theta_{t}, and corresponding geometric
% means.
alpha = vAlpha(t-lower_lim+1);
theta = vTheta(t-lower_lim+1);
GM_fx = vGM_fx(t-lower_lim+1);
GM_gx = vGM_gx(t-lower_lim+1);




end

