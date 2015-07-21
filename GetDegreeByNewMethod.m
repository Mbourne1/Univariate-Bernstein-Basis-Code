function [deg_calc,out_subresultants_unprocessed,out_subresultants_preprocessed,out_alphas,out_thetas,out_gm_fxs,out_gm_gxs] = ...
    GetDegreeByNewMethod(fx,gx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree of the AGCD of input polynomials fx and gx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs.

% fx : coefficients of polynomial f, expressed in Bernstein Basis gx :
% coefficients of polynomail g, expressed in Bernstein Basis

%                           Outputs.

% degree_calc - The calculated degree by various methods

% out_subresultants_unprocessed - All unprocesed subresultants S_{k} for k
% = 1,...,min(m,n)

% out_subresultants_preprocessed - All processed subresultants S_{k} for k
% = 1,...,min(m,n)

% out_alphas - All calculated optimal values of alpha from S_{k} for k =
% 1,...,min(m,n)

% out_thetas - All calculated optimal values of theta from S_{k} for k =
% 1,...,min(m,n)

% out_gm_fxs - All calculated geometric means of each C_{k}(f)

% out_gm_gxs - All calculated geometric means of each C_{k}(g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BOOL_PREPROC - (Boolean) It has been shown that the inclusion of
% preprocessors Geometric mean, scaling by alpha, change of independent
% variable, yield improved results.
%   1 :- Include Preprocessors. 0 :- Exclude Preprocessors.
global bool_preproc

% plotgraphs (bool)
%   0 - Don't plot graphs, just perform root finding operation. 1 - Plot
%   Graphs associated with calculating the GCD
global bool_plotgraphs

% bool_reordercols (bool)
%   0 - Leave columns of the Sylvester Matrix as standard partitions, where
% the first n-k+1 columns contain entries corresponding to the coefficients
% of f(y), and the last m-k+1 columns correspond to coefficients of g(y)
%   1 - Rearrange columns of the Sylvester subresultant matrices in
% accordance with Z Zeng - Computing Multiple Roots of inexact polynomials
% (page 889)
global bool_reordercols

% geometricMeanMethod used when calculating geometric means of the entries
% of a Sylvester matrix, a standard method would consider each entry in the
% matrix, but a new method, described in my internal report offers speed up
% due to the structured nature of the Sylvester matrix entries for the
% Sylvester matrix in the Bernstein basis.
%   0 - use MatLab Built in method for calculating geometric means 1 - use
%   my method of calculating geometric means
global geometricMeanMethod

% nominal value used when: Only one sylvester subresultant exists, ie k = 1
% and min(m,n) = 1. where m and n are the degrees of input polynomials f
% and g. Or when all subresultants Sk for k = 1,...,min(m,n) are all rank
% deficient or all of full rank.
global nominal_value

global min_delta_mag_rowsum

switch bool_plotgraphs
    case 1
        PlotNormalisedDiags = 1;
        PlotNormalisedRowSums = 1;
        PlotQRMaxMinRowSum = 1;
        PlotQRMaxMinDiagElem = 1;
        PlotAlpha = 1;
        PlotTheta = 1;
        PlotHeatMapSk = 0;
        PlotHeatMapQR = 0;
        PlotResiduals = 1;
    case 0
        PlotNormalisedDiags = 0;
        PlotNormalisedRowSums = 0;
        PlotQRMaxMinRowSum = 0;
        PlotQRMaxMinDiagElem = 0;
        PlotAlpha = 0;
        PlotTheta = 0;
        PlotHeatMapSk = 0;
        PlotHeatMapQR = 0;
        PlotResiduals = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get degree m of polynomial f
m = length(fx)-1;

% Get degree n of polynomial g
n = length(gx)-1;

% get minimum degree of f and g
min_mn = min(m,n);

% Initialise vectors to store all optimal alphas and theta, and each
% geometric mean for f and g in each S_{k} for k = 1,...,min(m,n)
alpha_vec    =   zeros(1,min_mn);
theta_vec    =   zeros(1,min_mn);

gm_fx_vec    =   zeros(1,min_mn);
gm_gx_vec    =   zeros(1,min_mn);

% Initialise vectors to store values calculated from each subresultant
% S_{k} for k = 1,...,min(m,n).

% Initialise a vector to store minimal residuals obtained by QR
% decomposition of each subresultant S_{k} for k=1,...,min(m,n)
minResQR_vec =   zeros(1,min_mn);

% Initialise a vector to store minimal residuals obtained by SVD of each 
% subresultant S_{k} for k=1,...,min(m,n)
minResSVD_vec =  zeros(1,min_mn);

% Initialise a vector to store max/min diagonal entry in the upper
% triangular matrix R1_{k} from the QR decomposition of S_{k} for 
% k=1,...,min(m,n)
ratio_maxmin_diag_vec = zeros(1,min_mn);

% Initialise a vector to store max row sum / min row sum of the rows of
% R1_{k} from the QR decomposition of S_{k} for k=1,...,min(m,n)
ratio_maxmin_rowsum_vec  = zeros(1,min_mn);

% Stores Data from QR decomposition.
Data_RowNorm = [];
Data_DiagNorm = [];

% Initialise arrays to store each subresultant matrix, both processed and
% unprocessed.
Sylvester_array_unproc = cell(min_mn,1);
Sylvester_array_preproc = cell(min_mn,1);
Cf = cell(min_mn,1);
Cg = cell(min_mn,1);


% For each subresultant $$S_{k}$$
for k = 1:1:min_mn
    
    switch bool_preproc
        case 1 % Include Preprocessors
            % if the first subresultant, then build from scratch
            if (k==1)
                % Include Preprocessor.
                C_f_unproc = BuildToeplitz(fx,1,n,k);
                C_g_unproc = BuildToeplitz(gx,1,m,k);
                
                
                % Get Geometric means.
                switch geometricMeanMethod
                    case 0 % use matlab method
                        % Get Geometric mean matlab method
                        gm_fx_vec(k) = geomean(abs(C_f_unproc(C_f_unproc~=0)));
                        gm_gx_vec(k) = geomean(abs(C_g_unproc(C_g_unproc~=0)));
                        
                    case 1 % use my method
                        % Get Geometric means my method
                        gm_fx_vec(k) = GeometricMean(abs(fx),n,k);
                        gm_gx_vec(k) = GeometricMean(abs(gx),m,k);
                    case 2 
                        % Do not calculate GeometricMean
                        gm_fx_vec(k) = 1;
                        gm_gx_vec(k) = 1;
                    case 3
                        gm_fx_vec(k) = mean(fx);
                        gm_gx_vec(k) = mean(gx);
                end
                
                % Divide normalised polynomials in Bernstein basis by
                % geometric means.
                fx_n = fx./gm_fx_vec(k);
                gx_n = gx./gm_gx_vec(k);
                
                % Divide the Sylvester Matrix partitions by Geometric mean.
                C_f_unproc = C_f_unproc ./ gm_fx_vec(k);
                C_g_unproc = C_g_unproc ./ gm_gx_vec(k);
                
                
            else % all subsequent subresultants, build from new method of
                %'build up' A_{k}S_{k}B_{k}
                
                % Get unprocessed partitions
                C_f_unproc = BuildToeplitz2(m,n,k-1,C_f_unproc);
                C_g_unproc = BuildToeplitz2(n,m,k-1,C_g_unproc);
                
                
                % Calculate Geometric mean. Note: my method depends on
                % using from scratch method. So set geometric mean Method
                % to Matlab Method.
                %geometricMeanMethod = 0;
                switch geometricMeanMethod
                    case 0
                        % Get Geometric means
                        gm_fx_vec(k) = geomean(abs(C_f_unproc(C_f_unproc~=0)));
                        gm_gx_vec(k) = geomean(abs(C_g_unproc(C_g_unproc~=0)));
                        
                    case 1
                        % Get Geometric means
                        gm_fx_vec(k) = GeometricMean(abs(fx),n,k);
                        gm_gx_vec(k) = GeometricMean(abs(gx),n,k);
                        
                end
                
                % Divide by Geometric means
                C_f_unproc = C_f_unproc ./ gm_fx_vec(k);
                C_g_unproc = C_g_unproc ./ gm_gx_vec(k);
                
                % Divide normalised polynomials in Bernstein basis by
                % geometric means.
                fx_n = fx./gm_fx_vec(k);
                gx_n = gx./gm_gx_vec(k);
            end
            
            % Build subresultant S_{k}, and add to array of Sk
            Sylvester_array_unproc{k} = [C_f_unproc C_g_unproc];
            
            % For each coefficient ai of F, obtain the max and min such
            % that F_max = [max a0, max a1,...] and similarly for F_min,
            % G_max, G_min
            
            [F_max,F_min,G_max,G_min] = GetMaxMin(Sylvester_array_unproc{k},m,n,k);
            
            % Calculate the optimal value of alpha and theta for the kth
            % subresultant matrix.
            [alpha_vec(k),theta_vec(k)] = OptimalAlphaTheta(F_max,F_min,G_max,G_min);
            
            % If linprog fails to find alpha and theta, and we have
            % assigned value 1. Best approximation is previous value of
            % alpha and theta when it exists.
            try
                if (alpha_vec(k) == 0 || alpha_vec(k) == 1)
                    alpha_vec(k) = alpha_vec(k-1);
                end
                
                if (theta_vec(k) == 1 || theta_vec(k) == 0)
                    theta_vec(k) = theta_vec(k-1);
                end
                
            catch
                alpha_vec(k) = 1;
                theta_vec(k) = 1;
            end
            
            
        case 0
            % Exclude preprocessors
            
            % Dont normalise by geometric mean.
            fx_n = fx;
            gx_n = gx;
            % Don't obtain optimal alpha and theta.
            alpha_vec(k) = 1;
            theta_vec(k) = 1;
            gm_fx_vec(k) = 1;
            gm_gx_vec(k) = 1;
    end
    
    %
    % Calculate the coefficients of the modified Bernstein basis
    % polynomials F2 and G2. Multiply G2 by alpha.
    fw_n =  fx_n*(theta_vec(k).^(0:1:m)) ;%\bar{f}(\theta,\omega)
    gw_n =  gx_n*(theta_vec(k).^(0:1:n)) ;%\bar{g}(\theta,\omega)
    
    % Construct the kth subresultant matrix for the optimal values of alpha
    % and theta.
    
    if (k==1)
        % first subresultant must be constructed in the conventional sense.
        Cf{k} = BuildToeplitz(fx_n,theta_vec(k),n,k);
        Cg{k} = BuildToeplitz(gx_n,theta_vec(k),m,k);
        
        % Build subresultant
        Sk = [Cf{k} alpha_vec(k).*Cg{k}]   ;
        
    else
        % subsequent subresultants can be built using my build up method.
        Cf{k} = BuildToeplitz2(m,n,k-1,Cf{k-1});
        Cg{k} = BuildToeplitz2(n,m,k-1,Cg{k-1});
        Sk = [Cf{k} Cg{k}];
        
    end
    
    
    if bool_reordercols == 1 % roots example
        
        Sk_experiment = zeros(size(Sk));
        % for each column in sk assign a position
        for i = 1:1:n-k+1
            Sk_experiment(:,2*i) = Sk(:,i);
        end
        
        % for each column in Sk second partition
        for i = n-k+2:1:m+n-(2*k)+2
            Sk_experiment(:,(i-(n-k+1))*2-1) = Sk(:,i);
        end
        
        Sk = Sk_experiment;
    end
    
    Sylvester_array_preproc{k} = Sk;
    

    % Using QR Decomposition of the sylvester matrix
    [~,R] = qr(Sk);
    
    % Take absolute values.
    R = abs(R);
    
    % Get number of rows in R1
    [R1_rows,~] = size(diag(R));
    
    % Obtain R1, the top square of the R matrix
    R1 = R(1:R1_rows,1:R1_rows);
    
    % Plot a heat map of the values of the R1 matrix
    switch PlotHeatMapQR
        case 1
            figure(18);
            colormap('hot');
            imagesc(log10(R1));
            caxis([-20,2])
            colorbar;
            %saveas(gcf,sprintf('Outputs/image(%i).jpg',k),'jpg')
    end
    
    % Get the diagonal elements of rows of R1
    R1_diag = abs(diag(R));
    
    % Get Norms of each row in matrix R1.
    R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);
    
    % Get ONLY the diagonal elements and normalise them.
    R1_DiagNorm = diag(R1)./norm(diag(R1));
    
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
    ratio_maxmin_diag_vec(k) = max(diag(R1))./min(diag(R1));
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
    
    % For each subresultant Sk - Remove each column c_{k,i} in turn, obtain
    % the residuals c_{k}-A_{k}.
    
    % residualQR_vector :- stores all residuals for 1,...,m+n-k+2 for a
    % given subresultant
    
    % minResQR_vector :- stores only one residual (the min) for each
    % subresultant k=1,...min(m,n)
    
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
    
    % Note - We could return the optimal column for each subresultant S_{k}
    % here, then when degree t of GCD is obtained, use this for optimal
    % column in calculating coprime polynomials u and v. At the moment,
    % optimal column is in a seperate file, called from o1.m
    
end

% % Edit 06/05/2015 Check to see if all subresultants are rank deficient in
% which case the degree of the GCD is min(m,n)


delta_mag_rowsum = max(abs(diff(log10(ratio_maxmin_rowsum_vec))));
%plot(log10(ratio_maxmin_rowsum_vec),'-s')

% % Edit 06/05/2015
% Check to see if only one subresultant exists, ie if m or n is equal
% to one

if min_mn == 1
    fprintf('Degree of GCD is either one or zero\n')

    max_r = max(abs(log10(R1_RowNorm)));
    min_r = min(abs(log10(R1_RowNorm)));
    
  
    if max_r./min_r > nominal_value
        deg_calc = 1;
        out_subresultants_unprocessed = Sylvester_array_unproc;
        out_subresultants_preprocessed = Sylvester_array_preproc;
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
    
    % get row sums
    %         figure()
    %         hold on
    %         plot(R1_RowNorm)
    %         hold off
    
elseif abs(delta_mag_rowsum) < min_delta_mag_rowsum
    % % Edit 06/05/2015
    % Check to see if all subresultants are rank deficient in which case
    % the degree of the GCD is min(m,n)
    
    %x = sum(diff(abs(log10(ratio_maxmin_rowsum_vector))))
    
    %sum(diff(abs(log10(ratio_maxmin_diag_vector))))
    
    fprintf('\n')
    fprintf('All subresultants appear to be either rank defficient or of full rank \n')
    fprintf('degree of GCD is either equal to min(m,n) or zero \n')
    fprintf('\n')
    
    deg_calc = min(m,n);
    fprintf('--------------------------------------------------------------------------- \n')
    fprintf('Degree By "From Scratch" Method.\n\n')
    fprintf('Calculated Degree of AGCD: %i \n', deg_calc)
    fprintf('--------------------------------------------------------------------------- \n')
else
    
    % Check to see if all subresultants are full rank, in which case the
    % degree of the GCD is 0
    
    
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
    % Plot Graph of ratio of max min elements.
    switch PlotQRMaxMinDiagElem
        case 1
            % Plot Graph of ratio of max : min element of the diagonal
            % elements of R1 from the QR decompositions.
            figure(11)
            x = 1:min_mn;
            plot(x,log10(ratio_maxmin_diag_vec),'red-s');
            axis([1,min_mn,0,inf])
            hold on
            legend('max/min diag element of subresultant S_{k}');
            title('max:min diagonal elements of R1 from The QR decomposition of S_{k} \newline (New Method)');
            hold off
    end
    
    % Plot Graph of ratio of max : min row sum in R1 from the QR
    % decompositions.
    switch PlotQRMaxMinRowSum
        case 1
            % Plot Graph of ratio of max : min row sum in R1 from the QR
            % decompositions.
            figure(12)
            x = 1:min_mn;
            plot(x,log10(ratio_maxmin_rowsum_vec),'red-s');
            hold on
            axis([1,min_mn,0,inf])
            legend('max/min Row Sum of Rows in R1 from the QR decomposition of S_{k}');
            title('Max:Min Row sum of Rows in R1 from the QR Decomposition of S_{k} \newline (New Method)');
            hold off
    end
    
    % Plot graph of norms of each row (N) from the qr decompostion of each
    % S_{k}
    switch PlotNormalisedRowSums
        case 1
            figure(13)
            plot(Data_RowNorm(:,1),(log10(Data_RowNorm(:,2))),'*')
            axis([0.9,min_mn,-inf,+inf])
            xlabel('k')
            ylabel('Normalised Row Sums of R1 in S_{k}')
            title(['Normalised Row Sums of R1 fom the QR decomposition of each subresultant S_{k} \newline '...
                'm = ' int2str(m) ', n = ' int2str(n) ' (New Method)']);
            hold off
    end
    
    switch PlotNormalisedDiags
        case 1
            figure(14)
            plot(Data_DiagNorm(:,1),(log10(Data_DiagNorm(:,2))),'*')
            axis([0.9,min_mn,-inf,+inf])
            xlabel('k')
            ylabel('Normalised Diagonals of R1 in S_{k}')
            title(['Normalised Diagonals in R1 matrix from the QR decomposition of each subresultant S_{k} \newline '...
                'm = ' int2str(m) ', n = ' int2str(n) ' \newline (New Method)']);
            hold off
    end
    
    
    % Plot graph of values of alpha for each subresultant
    switch PlotAlpha
        case 1
            figure(15)
            plot(1:1:length(alpha_vec),log10(alpha_vec),'-s')
            hold on
            xlabel('k')
            ylabel('log_{10} \alpha')
            title('Optimal values of \alpha for each subresultant S_{k} \newline (New Method)')
            hold off
    end
    
    % Plot graph of values of theta for each subresultant
    
    switch PlotTheta
        case 1
            figure(16)
            plot(1:1:length(theta_vec),log10(theta_vec),'-s')
            hold on
            xlabel('k')
            ylabel('log_{10} \theta')
            title('Optimal values of \theta for each subresultant S_{k} \newline (New Method)')
            hold off
            
    end
    
    
    % Plot graph of Residuals by QR and SVD, if using the residual method
    % to calculate the degree of the GCD.
    switch PlotResiduals
        case 1
            x = 1:1:k;
            figure(17)
            plot(x,log10(minResQR_vec),'red-o','DisplayName','Residuals by QR');
            hold on
            plot(x,log10(minResSVD_vec),'blue-s','DisplayName','Residauls by SVD');
            axis([1,min_mn,-inf,+inf])
            ylabel('log_{10} Residual')
            xlabel('k')
            title('Residual obtained by removing optimal column from S_{k} \newline (New Method)');
            legend(gca,'show')
            hold off
    end

 % Outputs
    
    % Output just corresponding to calculated value of the degree. Output
    % subresultant S_{t}, alpha_{t}, theta_{t}, and corresponding geometric
    % means.
    out_subresultant_unprocessed = cell2mat(Sylvester_array_unproc(deg_calc));
    out_subresultant_preprocessed = cell2mat(Sylvester_array_preproc(deg_calc));
    out_alpha = alpha_vec(deg_calc);
    out_theta = theta_vec(deg_calc);
    out_gm_fx = gm_fx_vec(deg_calc);
    out_gm_gx = gm_gx_vec(deg_calc);
    
    % Output all subresultants, all optimal alphas, all optimal thetas and
    % all geometric means for each subresultant S_{k} where k =
    % 1,...,min(m,n)
    out_subresultants_unprocessed = Sylvester_array_unproc;
    out_subresultants_preprocessed = Sylvester_array_preproc;
    out_alphas = alpha_vec;
    out_thetas = theta_vec;
    out_gm_fxs = gm_fx_vec;
    out_gm_gxs = gm_gx_vec;


end



