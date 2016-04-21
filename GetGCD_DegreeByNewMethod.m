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

% Global Variables
global BOOL_ALPHA_THETA
global PLOT_GRAPHS
global GEOMETRIC_MEAN_METHOD
global THRESHOLD



% Get degree m of polynomial f
m = GetDegree(fx);

% Get degree n of polynomial g
n = GetDegree(gx);

% get minimum degree of f and g
min_mn = min(m,n);

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(1,min_mn);
vTheta    =   zeros(1,min_mn);

vGM_fx    =   zeros(1,min_mn);
vGM_gx    =   zeros(1,min_mn);

% Initialise vectors to store values calculated from each subresultant
% S_{k} for k = 1,...,min(m,n).

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
minResQR_vec =   zeros(1,min_mn);

% Initialise a vector to store minimal residuals obtained by SVD of each
% subresultant S_{k} for k=1,...,min(m,n)
minResSVD_vec =  zeros(1,min_mn);

vMaxDiagR1= zeros(1,min_mn);
vMinDiagR1= zeros(1,min_mn);
vMaxRowNormR1= zeros(1,min_mn);
vMinRowNormR1 = zeros(1,min_mn);

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
    
    [vGM_fx(k), vGM_gx(k), vAlpha(k), vTheta(k)] = Preprocess(fx,gx,k);
    
    % Divide the Sylvester Matrix partitions by Geometric mean.
    DT1Q1_unproc = DT1Q1_unproc ./ vGM_fx(k);
    DT2Q2_unproc = DT2Q2_unproc ./ vGM_gx(k);
    
    
    % Divide normalised polynomials in Bernstein basis by
    % geometric means.
    fx_n = fx./vGM_fx(k);
    gx_n = gx./vGM_gx(k);
    
    %
    % Calculate the coefficients of the modified Bernstein basis
    % polynomials F2 and G2. Multiply G2 by alpha.
    
    fw_n = GetWithThetas(fx_n,vTheta(k));
    gw_n = GetWithThetas(gx_n,vTheta(k));
    
    % Construct the kth subresultant matrix for the optimal values of alpha
    % and theta.
    
    if (k==1)
        % first subresultant must be constructed in the conventional sense.
        DT1Q1 = BuildDT1Q1(fw_n,n-k);
        DT2Q2 = BuildDT1Q1(gw_n,m-k);
        
        % Build subresultant
        DTQ = [DT1Q1 vAlpha(k).*DT2Q2]   ;
        
    else
        % subsequent subresultants can be built using my build up method.
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
    
    % For each subresultant Sk - Remove each column c_{k,i} in turn, obtain
    % the residuals c_{k}-A_{k}.
    
    % residualQR_vector :- stores all residuals for 1,...,m+n-k+2 for a
    % given subresultant
    
    % minResQR_vector :- stores only one residual (the min) for each
    % subresultant k=1,...min(m,n)
    
    vMinimumResidual(k) = GetMinimalDistance(DTQ);
    
    
end

vRatio_MaxMin_Diagonals_R = vMaxDiagR1 ./ vMinDiagR1;
vRatio_MaxMin_RowNorm_R = sanitize(vRatio_MaxMin_Diagonals_R);

vRatio_MaxMin_RowNorm_R = vMaxRowNormR1 ./ vMinRowNormR1;
vRatio_MaxMin_RowNorm_R = sanitize(vRatio_MaxMin_RowNorm_R);



% Get the change in the ratios from one subresultant to the next.
vDelta_MaxMin_RowNorm_R = abs(diff(log10(vRatio_MaxMin_RowNorm_R)));

% Get the change in the ratios of diagonal elements from one subresultant
% to the next.
vDelta_MaxMin_Diag_R = abs(diff(log10(vRatio_MaxMin_Diagonals_R)));

% Get the maximum change in rowsum ratio and its index
[max_delta_mag_rowsum,indexMaxChange_RowNorm] = max(vDelta_MaxMin_RowNorm_R);

% Get the maximum change in diag ratio and its index
[max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = max(vDelta_MaxMin_Diag_R);


% Check to see if only one subresultant exists, ie if m or n is equal
% to one

if min_mn == 1
    t = Get_Rank_One_Subresultant(DTQ);
    alpha = vAlpha(1);
    theta = vTheta(1);
    GM_fx = vGM_fx(1);
    GM_gx = vGM_gx(1);
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
        alpha = 0;
        theta = 0;
        GM_fx = 0;
        GM_gx = 0;
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

%% Graph Plotting

switch PLOT_GRAPHS
    case 'y'
        % Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
        figure('name',figure_name)
        x = 1:min_mn;
        plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
        hold on
        axis([1,min_mn,0,inf])
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
        
        
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure_name = sprintf('%s : Max:min Row Norms',mfilename);
        figure('name',figure_name)
        x = 1:min_mn;
        plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
        hold on
        axis([1,min_mn,0,inf])
        legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
        hold off
        
        
        % Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
        figure_name = sprintf('%s : RowNorm',mfilename);
        figure('name',figure_name)
        plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
        axis([0.9,min_mn,-inf,+inf])
        xlabel('k')
        ylabel('Normalised Row Sums of R1 in S_{k}')
        title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
        
        
        % %
        figure_name = sprintf('%s : Diagonals',mfilename);
        figure('name',figure_name)
        plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
        axis([0.9,min_mn,-inf,+inf])
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


% Outputs

% Output just corresponding to calculated value of the degree. Output
% subresultant S_{t}, alpha_{t}, theta_{t}, and corresponding geometric
% means.
alpha = vAlpha(t);
theta = vTheta(t);
GM_fx = vGM_fx(t);
GM_gx = vGM_gx(t);




end



