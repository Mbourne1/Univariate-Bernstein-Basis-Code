function [t,alpha, theta,GM_fx,GM_gx] = ...
    GetGCD_Degree(fx,gx)
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
% t : Degree of GCD of f(x) and g(x) 
%
% alpha : optimal value of alpha
%
% theta : optimal value of theta
% 
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
m = GetDegree(fx);

% Get degree n of polynomial g(x)
n = GetDegree(gx);

fprintf('Begin : Get Degree m = % i , n = % i \n\n' ,m,n)

% get minimum degree of f(x) and g(x)
min_mn = min(m,n);

%% Initialisation stage

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
vAlpha    =   zeros(min_mn,1);
vTheta    =   zeros(min_mn,1);
vGM_fx    =   zeros(min_mn,1);
vGM_gx    =   zeros(min_mn,1);

vMax_Diag_R1 = zeros(min_mn,1);
vMin_Diag_R1 = zeros(min_mn,1);
vMax_R1_RowNorm = zeros(min_mn,1);
vMin_R1_RowNorm = zeros(min_mn,1);

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
vMinimumResidual             =   zeros(min_mn,1);

% Initialise a vector to store max/min diagonal entry in the upper
% triangular matrix R1_{k} from the QR decomposition of S_{k} for
% k=1,...,min(m,n)
vRatio_MaxMin_Diagonals_R    =   zeros(min_mn,1);


vMinimumSingularValues = zeros(min_mn,1);

% Stores Data from QR decomposition.
Data_RowNorm    = [];
Data_DiagNorm   = [];


%% Loop

% For each subresultant $S_{k}$
for k = 1:1:min_mn
    
    
    [vGM_fx(k), vGM_gx(k),vAlpha(k),vTheta(k)] = Preprocess(fx,gx,k);
    
    
    % Divide f(x) and g(x) by geometric means
    fx_n = fx./vGM_fx(k);
    gx_n = gx./vGM_gx(k);
    
    % Construct the kth subresultant matrix S_{k}(f(\theta),g(\theta))
    fw = GetWithThetas(fx_n,vTheta(k));
    gw = GetWithThetas(gx_n,vTheta(k));
    
    % Build the k-th subresultant
    Sk = BuildSubresultant(fw,vAlpha(k).*gw,k);
    
    % Get the minimal singular value from S_{k}
    vMinimumSingularValues(k) = min(svd(Sk));
    
    % Get the matrix R1 from the QR Decomposition of S
    R1 = GetR1(Sk);
    
    % Get Norms of each row in the matrix R1
    vR1_RowNorms = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get Norms of diagonal eleements of R1
    vR1_DiagNorm = diag(R1)./norm(diag(R1));
    
    % Insert Row Norm data into a matrix
    Data_RowNorm = AddToResults(Data_RowNorm,vR1_RowNorms,k);
    
    % Insert diagonals data into a matrix
    Data_DiagNorm = AddToResults(Data_DiagNorm,vR1_DiagNorm,k);
    
    % Get maximum and minimum row diagonals of R1
    vMax_Diag_R1(k) = max(diag(R1));
    vMin_Diag_R1(k) = min(diag(R1));
       
    % Get maximum and minimum row norms of rows of R1.
    vMax_R1_RowNorm(k) = max(vR1_RowNorms);
    vMin_R1_RowNorm(k) = min(vR1_RowNorms);
    
    % Get the minimal residuals for each subresultant S_{k}.
    vMinimumResidual(k) = GetMinimalDistance(Sk);
    
    
end

% End of Loop

% Get ratio of max:min Row Norms of R1
vRatio_MaxMin_RowNorm_R = vMax_R1_RowNorm./ vMin_R1_RowNorm;
vRatio_MaxMin_RowNorm_R = sanitize(vRatio_MaxMin_RowNorm_R);

% Get ratio of max:min diagonals of R1
vRatio_MaxMin_Diagonals_R = vMax_Diag_R1 ./ vMin_Diag_R1;
vRatio_MaxMin_Diagonals_R = sanitize(vRatio_MaxMin_Diagonals_R);



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

%% Analysis of Minimum Singular values

%%

% If only one subresultant exists, get GCD by alternative method
if min_mn == 1
    t = Get_Rank_One_Subresultant(Sk);
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
figure('name','GetDegree - MaxMin - Row Diags')
x = 1:min_mn;
plot(x,log10(vRatio_MaxMin_Diagonals_R),'red-s');
hold on
axis([1,min_mn,0,inf])
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off


% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
figure('name','GetDegree - MaxMin - Row Norms')
x = 1:min_mn;
plot(x,log10(vRatio_MaxMin_RowNorm_R),'red-s');
hold on
axis([1,min_mn,0,inf])
legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');
title('Max:Min Row Norms of Rows in R1 from the QR Decomposition of S_{k} (Original)');
hold off


% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}

figure('name','GetDegree - RowNorm')
plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
axis([0.9,min_mn,-inf,+inf])
xlabel('k')
ylabel('Normalised Row Sums of R1 in S_{k}')
title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
    'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
hold off


figure('name','GetDegree - DiagNorm')
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


%%
% Outputs

% Output all subresultants, all optimal alphas, all optimal thetas and all
% geometric means for each subresultant S_{k} where k = 1,...,min(m,n)

alpha                   = vAlpha(t);
theta                   = vTheta(t);
GM_fx                   = vGM_fx(t);
GM_gx                   = vGM_gx(t);
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

