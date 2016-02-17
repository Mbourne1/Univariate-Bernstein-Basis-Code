function [deg_calc,out_subresultants_unprocessed,out_subresultants_preprocessed,out_alphas,out_thetas,out_gm_fxs,out_gm_gxs] = ...
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

global bool_preproc
global plot_graphs
global nominal_value
global min_delta_mag_rowsum


switch plot_graphs
    case 'y'
        PlotNormalisedDiagonals = 1;
        PlotNormalisedRowSums = 1;
        PlotQRMaxMinRowSum = 1;
        PlotQRMaxMinDiagElem = 1;
        PlotAlpha = 1;
        PlotTheta = 1;
        PlotResiduals = 1;
    case 'n'
        PlotNormalisedDiagonals = 0;
        PlotNormalisedRowSums = 0;
        PlotQRMaxMinRowSum = 0;
        PlotQRMaxMinDiagElem = 0;
        PlotAlpha = 0;
        PlotTheta = 0;
        PlotResiduals = 0;
    otherwise
        error('bool_plotgraphs is either y or n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get degree m of polynomial f(x)
[r,~] = size(fx);
m = r - 1;

% Get degree n of polynomial g(x)
[r,~] = size(gx);
n = r - 1;

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
minResQR_vec             =   zeros(min_mn,1);

% Initialise a vector to store minimal residuals obtained by SVD of each
% subresultant S_{k} for k=1,...,min(m,n)
minResSVD_vec            =   zeros(min_mn,1);

% Initialise a vector to store max/min diagonal entry in the upper
% triangular matrix R1_{k} from the QR decomposition of S_{k} for
% k=1,...,min(m,n)
ratio_maxmin_diag_vec    =   zeros(min_mn,1);

% Initialise a vector to store max row sum / min row sum of the rows of
% R1_{k} from the QR decomposition of S_{k} for k=1,...,min(m,n)
ratio_maxmin_rowsum_vec  =   zeros(min_mn,1);

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
    
    % Whether or not applying Preprocessors
    
    % Get Unprocessed partitions (Including Geometric Mean)
    C_f_unproc = BuildT1(fx,1,n,k);
    C_g_unproc = BuildT1(gx,1,m,k);
    
    switch bool_preproc
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
    
    % Construct the kth subresultant matrix for the optimal
    % values of alpha and theta.
    Sk  =   BuildSubresultant(fx_n,gx_n,k,alpha_vec(k),theta_vec(k));
    
    % Add Sk to the array of preprocessed Sk
    Sylvester_array_preprocessed{k} = Sk;
    
    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(Sk);
    
    % Take absolute values.
    R = abs(R);
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1 the top square of the R matrix.
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Get Norms of each row in the matrix R1
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
    
    % Scatter Plot Data
    ks = k.*ones(size(R1_RowNorm));
    ns = 1:1:size(R1_RowNorm,1);
    
    % Form a triple of [ks, the value of QR_RowNorm, and the index of the value of
    % the row of R1 corresponding to QR_RowNorm].
    % EG.
    %  [1   0.015  1
    %   1   0.156  2
    %   2 ...]
    X = [ks R1_RowNorm ns'];
    Data_RowNorm = [Data_RowNorm; X];
    
    X2 = [ks R1_DiagNorm ns'];
    Data_DiagNorm = [Data_DiagNorm;X2];
    
    
    
    % Get max:min diag elem of S_{k}
    ratio_maxmin_diag_vec(k) = max(diag(R1))./min(diag(R1));
    
    % Get max:min rownorm r_{i}/r_{j} of S_{k}
    ratio_maxmin_rowsum_vec(k) = max(R1_RowNorm)./min(R1_RowNorm);
    
    % Edit - 23/02/2015 - remove any infinite values
    for i = length(ratio_maxmin_diag_vec):-1:2
        if isinf(ratio_maxmin_diag_vec(i-1)) ||  ratio_maxmin_diag_vec(i-1) == 0
            ratio_maxmin_diag_vec(i-1) = ratio_maxmin_diag_vec(i);
        end
    end
    % Edit 23/02/2015 - remove any infinite values
    for i = length(ratio_maxmin_rowsum_vec):-1:2
        if isinf(ratio_maxmin_rowsum_vec(i-1)) || ratio_maxmin_rowsum_vec(i-1) ==0
            ratio_maxmin_rowsum_vec(i-1) = ratio_maxmin_rowsum_vec(i);
        end
    end
    
    %
    
    % For each subresultant Sk - Remove each column c_{k,i} in
    % turn, obtain the residuals c_{k}-A_{k}.
    
    % residualQR_vector :- stores all residuals for 1,...,m+n-k+2
    % for a given subresultant
    
    % minResQR_vector :- stores only one residual (the min) for
    % each subresultant k=1,...min(m,n)
    
    residualQR_vector = zeros(1,m+n-2*k+2);
    residualSVD_vector = zeros(1,m+n-2*k+2);
    for i = 1:1:m+n-2*k+2
        [Ak,ck]   = RemoveSubresultantColumn(Sk,i);
        residualQR_vector(i)  = CalculateResidualQR(ck,Ak);
        residualSVD_vector(i) = CalculateResidualSVD(ck,Ak);
    end
    
    %get minimal residual for this subresultant
    minResQR_vec(k) = min(residualQR_vector);
    minResSVD_vec(k) = min(residualSVD_vector);
    
    % Note - We could return the optimal column for each
    % subresultant S_{k} here, then when degree t of GCD is
    % obtained, use this for optimal column in calculating coprime
    % polynomials u and v. At the moment, optimal column is in a
    % seperate file, called from o1.m
    
    
    
end

% Get the change in the ratios from one subresultant to the next.
delta_mag_rowsum = abs(diff(log10(ratio_maxmin_rowsum_vec)));

% Get the change in the ratios of diagonal elements from one subresultant
% to the next.
delta_mag_maxmin_diag = abs(diff(log10(ratio_maxmin_diag_vec)));

% Get the maximum change in rowsum ratio and its index
[max_delta_mag_rowsum,index] = max(delta_mag_rowsum);

% Get the maximum change in diag ratio and its index
[max_delta_mag_maxmin_diag, index2] = max(delta_mag_maxmin_diag);

% Get the second largest value
[second_lrgst_delta_mag_rowsum] = max(delta_mag_rowsum(delta_mag_rowsum<max(delta_mag_rowsum)));

% Get the second largest value
[second_lrgs_delta_mag_maxmin_diag] = max(delta_mag_maxmin_diag(delta_mag_maxmin_diag<max(delta_mag_maxmin_diag)));

switch plot_graphs
    case 'y'
        % Plot for report
        figure('name','Get Degree - Maxmin - Rowsums')
        hold on
        plot(log10(ratio_maxmin_rowsum_vec),'-s')
        xlabel('k : index of subresultant S')
        ylabel('log_{10} ratio max:min row sum r_{i} in R1')
        title('Plotting max:min row sum of the rows r_{i} of R1 from the QR decomposition of each subresultants S_{k}')
        hold off
        
        % Plot for report
        figure('name','Get Degree - Maxmin - Diags')
        hold on
        plot(log10(ratio_maxmin_diag_vec),'-s')
        xlabel('k : index of subresultant S_{k}');
        ylabel('log_{10} ratio max:min diagon entries of R1')
        title('Plotting max:min diagonal entries of R1 from the QR decomposition of each subresultants S_{k}')
        hold off
    case 'n'
    otherwise
        error('error:')
end

% Check to see if only one subresultant exists, ie if m or n is equal
% to one
if min_mn == 1
    fprintf('############## Exception ###################################\n\n')
    fprintf('min(m,n) = 1 \n')
    fprintf('Only One Sylvester Subresultant Exists \n')
    fprintf('Degree of GCD is either one or zero \n')
    
    figure(997)
    hold on
    plot(log10(diag(R1)),'-s')
    xlabel('i: index of the ith diagonal')
    ylabel('log10 Diagonals of R1 ')
    title('Diagonal values in R1 from QR Decomposition of S_{1}(From Scratch)')
    hold off
    
    % Get the changes in R1
    delta_log_diags_R1 = abs(diff(log10(diag(R1))));
    
    % Get the maximum change
    max_delta_log_diags_R1 = max(delta_log_diags_R1);
    
    % Plot a graph of all the Row norms for the Subresultant S_{1}
    figure(998)
    hold on
    plot(log10(R1_RowNorm),'-s')
    title('Norms of the rows of the R matrix from the QR decompositino of Subresultant S_{1}')
    xlabel('i: index of the ith row')
    ylabel('log10 row norms of R1')
    hold off
    
    % Get the changes in row sums
    delta_log_rowsum_R1 = abs(diff(log10(R1_RowNorm)));
    % Get the maximum change
    max_delta_log_rowsum_R1 = max(delta_log_rowsum_R1);
    
    fprintf('Max change in row sums for each r_{i} is given by : %i \n', max_delta_log_rowsum_R1)
    fprintf('Current Nominal Value : %i \n',nominal_value)
    
    
    
    fprintf('############## Exception ###################################\n')
    
    
    % if max/min is greater than nominal value, then we suppose that the
    % subresultant S_{1} is rank deficient, so degree of GCD is one.
    
    if max_delta_log_rowsum_R1 > nominal_value
        % The maximum change in row sum (delta) is significant
        % Subresultant is rank deficient
        % Set degree of GCD = 1.
        
        deg_calc = 1;
        
        out_subresultants_unprocessed = Sylvester_array_unproc;
        out_subresultants_preprocessed = Sylvester_array_preprocessed;
        out_alphas = alpha_vec;
        out_thetas = theta_vec;
        out_gm_fxs = gm_fx_vec;
        out_gm_gxs = gm_gx_vec;
        
        fprintf('--------------------------------------------------------------------------- \n')
        fprintf('Degree By "From Scratch" Method.\n\n')
        fprintf('Calculated Degree of AGCD: %i \n', deg_calc)
        fprintf('--------------------------------------------------------------------------- \n')
        return
        
        
    else
        % The maximum change in row sum (delta) is insignificant
        % Subresultant is of full rank
        % Set degree of GCD = 0
        
        % Set Degree of GCD to zero
        deg_calc = 0;
        out_subresultants_unprocessed = 0;
        out_subresultants_preprocessed = 0;
        out_alphas = 0;
        out_thetas = 0;
        out_gm_fxs = 0;
        out_gm_gxs = 0;
        
        
        fprintf('--------------------------------------------------------------------------- \n')
        fprintf('Degree By "From Scratch" Method.\n\n')
        fprintf('Calculated Degree of AGCD: %i \n', deg_calc)
        fprintf('--------------------------------------------------------------------------- \n')
        
        return
    end
    
    
    % Set a condition for which we consider the maximum change in row sums to
    % significant or insignificant.
elseif abs(max_delta_mag_rowsum) < min_delta_mag_rowsum
    
    
    % Check to see if all subresultants are rank deficient in which case
    % the degree of the GCD is min(m,n)
    
    fprintf('\n')
    fprintf('############## Exception ###################################\n')
    fprintf('All subresultants appear to be either rank defficient or of full rank \n')
    fprintf('Degree of GCD is either equal to min(m,n) or zero \n')
    fprintf('\n')
    fprintf('Delta : %i \n', abs(delta_mag_rowsum));
    fprintf('nominal value (min_delta_mag_rowsum : %i \n',min_delta_mag_rowsum);
    fprintf('############## Exception ###################################\n')
    
    
    
    
    if min(log10(ratio_maxmin_rowsum_vec)) < 10
        
        %all subresultants are full rank
        deg_calc = 0;
        
        out_subresultants_unprocessed = 0;
        out_subresultants_preprocessed = 0;
        out_alphas = 0;
        out_thetas = 0;
        out_gm_fxs = 0;
        out_gm_gxs = 0;
        fprintf('--------------------------------------------------------------------------- \n')
        fprintf('Degree By "From Scratch" Method.\n\n')
        fprintf('Calculated Degree of AGCD: %i \n', deg_calc)
        fprintf('--------------------------------------------------------------------------- \n')
        return;
    else
        % all subresultants are rank deficient
        deg_calc = min(m,n);
    end
    
    
    fprintf('--------------------------------------------------------------------------- \n')
    fprintf('Degree By "From Scratch" Method.\n\n')
    fprintf('Calculated Degree of AGCD: %i \n', deg_calc)
    fprintf('--------------------------------------------------------------------------- \n')
else
    
    
    
    
    % Otherwise the degree of the GCD is somewhere between.
    
    [~,degree_calc_1] = max(abs(diff(log10(ratio_maxmin_rowsum_vec))));
    
    % Get degree by ratio of max/min element in N, where N is a vector of the
    % norms of each row of R1.
    [~,degree_calc_2] = max(abs(diff(log10(ratio_maxmin_diag_vec))));
    
    % Get degree by minimal residual obtained by removing each column from each
    % subresultant, giving Ak x = ck. When residual is 'close to zero'. Matrix is rank
    % deficient. (Note if zero, log10(0) = inf. so use fudge factor)
    A = log10(minResQR_vec);
    
    [~,degree_calc_3] = max(abs(diff(A)));
    
    % Take the mode of the three degree calculation methods
    deg_calc = mode([degree_calc_1,degree_calc_2,degree_calc_3]);
    fprintf('--------------------------------------------------------------------------- \n')
    fprintf('Degree By "From Scratch" Method.\n\n')
    fprintf('GCD Degree by max:min row sum in R1 from QR decomposition : %i \n',degree_calc_1);
    fprintf('GCD Degree by max:min diagonal elements in R1 from QR decomposition : %i \n',degree_calc_2)
    fprintf('GCD Degree by residual using QR Decomposition : %i \n',degree_calc_3)
    fprintf('--------------------------------------------------------------------------- \n')
    fprintf('Calculated Degree of AGCD by Standard Method: %i \n', deg_calc)
    fprintf('--------------------------------------------------------------------------- \n')
    
    
    
    
    
end

%% Graph Plotting
% Plot Graph of ratio of max min elements.
switch PlotQRMaxMinDiagElem
    case 1
        % Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
        figure('name','GetDegree - MaxMin - Row Diags')
        x = 1:min_mn;
        plot(x,log10(ratio_maxmin_diag_vec),'red-s');
        hold on
        axis([1,min_mn,0,inf])
        legend('Max:Min diag element of subresultant S_{k}');
        title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
        ylabel('log_{10} max:min diag element')
        hold off
end

% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
switch PlotQRMaxMinRowSum
    case 1
        % Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
        figure('name','GetDegree - MaxMin - Row Sums')
        x = 1:min_mn;
        plot(x,log10(ratio_maxmin_rowsum_vec),'red-s');
        hold on
        axis([1,min_mn,0,inf])
        legend('Max:Min Row Sum of Rows in R1 from the QR decomposition of S_{k}');
        title('Max:Min Row sum of Rows in R1 from the QR Decomposition of S_{k} (Original)');
        hold off
end

% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
switch PlotNormalisedRowSums
    case 1
        figure('name','GetDegree - RowNorm')
        plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
        axis([0.9,min_mn,-inf,+inf])
        xlabel('k')
        ylabel('Normalised Row Sums of R1 in S_{k}')
        title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
        
end
switch PlotNormalisedDiagonals
    case 1
        figure(4)
        plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
        axis([0.9,min_mn,-inf,+inf])
        xlabel('k')
        ylabel('Normalised Diagonals of R1 in S_{k}')
        title(['Normalised Diagonals in R1 matrix from the QR decomposition of each subresultant S_{k} \newline '...
            'm = ' int2str(m) ', n = ' int2str(n) '(Original)']);
        hold off
        
end

% Plot graph of values of alpha for each subresultant
switch PlotAlpha
    case 1
        figure('name','GetDegree - Alphas')
        plot(1:1:length(alpha_vec),log10(alpha_vec),'-s')
        hold on
        xlabel('k')
        ylabel('log_{10} \alpha')
        title('Optimal values of \alpha for each subresultant S_{k} (Original)')
        hold off
end

% Plot graph of values of theta for each subresultant
switch PlotTheta
    case 1
        figure('name','GetDegree - Thetas')
        plot(1:1:length(theta_vec),log10(theta_vec),'-s')
        hold on
        xlabel('k')
        ylabel('log_{10} \theta')
        title('Optimal values of \theta for each subresultant S_{k} (Original)')
        hold off
end

% Plot graph of Residuals by QR and SVD, if using the residual method to
% calculate the degree of the GCD.
switch PlotResiduals
    case 1
        x = 1:1:k;
        figure(7)
        plot(x,log10(minResQR_vec),'red-o','DisplayName','Residuals by QR');
        hold on
        plot(x,log10(minResSVD_vec),'blue-s','DisplayName','Residuals by SVD');
        axis([1,min_mn,-inf,+inf])
        ylabel('log_{10} Residual')
        xlabel('k')
        title('Residual obtained by removing optimal column from S_{k} (Original)');
        legend(gca,'show')
        hold off
end
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