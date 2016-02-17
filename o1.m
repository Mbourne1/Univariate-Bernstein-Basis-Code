function [f_output,g_output,d_output, u_output, v_output,alpha_opt,theta_opt] = ...
    o1(fx,gx)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
%                             Inputs:
%
%
% ex - (Int) Example Number
%
% emin - Signal to noise ratio (minimum)
%
% emax - Signal to noise ratio (maximum)
%
%
%                           Outputs:
%
% gcd_output -
%
% u_output -
%
% v_output -
%

%%
%                       GLOBAL VARIABLES

global bool_apf
global bool_q
global bool_sntln
global bool_preproc
global bool_denom_syl
global bool_SNTLN_Roots
global bool_APF_Roots
global plot_graphs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the degree m of polynomial f
m = length(fx) - 1;

% Get the degree n of polynomial g
n = length(gx) - 1;

%% Get degree of GCD by first method

% Method 1 - 'From Scratch' Build each subresultant from scratch
[t_byGetDegree, subres_unproc1, subres_preproc1, alphas_opt1, thetas_opt1, gms_fx1, gms_gx1] = ...
    GetDegree(fx,gx);

% %% Get Degree of GCD by second method
% % Method 2 - 'Build Up using S_{k} = A_{k}S_{k-1}B{k}' Build each
% % subresultant by matrix manipulation of the previous subresultant.
%
% switch bool_q
%     case 'n'
%         % If Q is excluded, we can not build up. since the form of the
%         % entries of S_k rely on the inclusion of Q.
%         t_byGetDegreeNew = t_byGetDegree;
%     case 'y'
%         [t_byGetDegreeNew, subres_unproc2, subres_preproc2, alphas_opt2, thetas_opt2, gms_fx2, gms_gx2] = ...
%             GetDegreeByNewMethod(fx,gx);
%     otherwise
%         error('bool_q must be set to either (y) or (n)')
% end
%
% if t_byGetDegreeNew == 0
%     fprintf('f(x) and g(x) appear to be coprime \n')
%     f_output = fx;
%     g_output = gx;
%     d_output = 1;
%     u_output = fx;
%     v_output = gx;
%     alpha_opt = 1;
%     theta_opt = 1;
%     return
% end


% Method 2 returns good results for obtaining the degree, but not good
% results for calculating u,v and GCD.


% switch bool_q
%     case 'n'
%         % if Q is excluded from the sylvester matrix, we were not able to
%         % use method 2, so take all values from method 1.
%         t = t_byGetDegree;
%         subresultant_unproc    = cell2mat(subres_unproc1(t));
%         subresultant_preproc   = cell2mat(subres_preproc1(t));
%         alpha_opt              = alphas_opt1(t);
%         theta_opt              = thetas_opt1(t);
%         gm_fx                  = gms_fx1(t);
%         gm_gx                  = gms_gx1(t);
%
%     case 'y'
%         % if Q was included in the sylvester matrix,
%
%         % we chose from method 1 or 2
%         method = 1;
%         switch method
%             case 1
%t = t_byGetDegreeNew;

t = t_byGetDegree;
subresultant_unproc    = cell2mat(subres_unproc1(t));
subresultant_preproc   = cell2mat(subres_preproc1(t));
alpha_opt              = alphas_opt1(t);
theta_opt              = thetas_opt1(t);
gm_fx                  = gms_fx1(t);
gm_gx                  = gms_gx1(t);


%             case 2
%                 t = t_byGetDegreeNew;
%                 subresultant_unproc    = cell2mat(subres_unproc2(t));
%                 subresultant_preproc   = cell2mat(subres_preproc2(t));
%                 alpha_opt              = alphas_opt2(t);
%                 theta_opt              = thetas_opt2(t);
%                 gm_fx                  = gms_fx2(t);
%                 gm_gx                  = gms_gx2(t);
%         end
%     otherwise
%         error('bool_q must be either y or n')
% end

% Always use the second method ('Build Up') Method to obtain the degree t
% of the GCD.

% If finding the GCD fails, set the degree of the GCD to be 1.
if isempty(t)
    t = 1;
end

% Normalise fx and gx by Geometric mean to obtain fx_n and gx_n Normalise
% by geometric mean obtained from subresultant S_{t}
switch bool_preproc
    case 'n' % Exclude preprocessors
        fx_normalised = fx;
        gx_normalised = gx;
        alpha_opt = 1;
        theta_opt = 1;
    case 'y' % include preprocessors
        fx_normalised = fx./gm_fx;
        gx_normalised = gx./gm_gx;
    otherwise
        error('bool_preproc must be either y or n')
end


%% Build a series of subresultants for analysis
% Sylvester_1 : Unprocessed
% Sylvester_2 : Preprocessed
% Sylvester_3 : With SNTLN

% Get Subresultant of unprocessed f(x) and g(x)
St_unproc = BuildSubresultant(fx,gx,1,1,1);

% Get Subresultant of preprocessed f(w) and g(w)
St_preproc = BuildSubresultant(fx_normalised,gx_normalised,1,alpha_opt,theta_opt);

%%
% Get the optimal column of the sylvester matrix to be removed. Where
% removal of the optimal column gives the minmal residual in (Ak x = ck)
[opt_col] = GetOptimalColumn(subresultant_preproc);

%% Perform SNTLN
% Apply / Don't Apply structured perturbations to Sylvester Matrix NOTE
% SNTLN SHOULD NOT BE USED WITH (METHOD 2 - 'BUILD UP'). SNTLN uses build
% from scratch method.

switch bool_sntln
    case 'y'
        % Obtain refined values of fx and gx, alpha, theta, and vector x
        % consisting of coefficients of u and v
        
        
        switch bool_SNTLN_Roots
            case 'StandardSNTLN' % use standard SNTLN
                [fx_n,gx_n,alpha_opt,theta_opt,X] = ...
                    SNTLN(fx_normalised,gx_normalised,alpha_opt,theta_opt,t,opt_col);
                
            case 'RootSpecificSNTLN' % use roots specific sntln
                [fx_n,gx_n,alpha_opt,theta_opt,X] = ...
                    SNTLN_Roots(fx_normalised,gx_normalised,alpha_opt,theta_opt,t,opt_col,gm_fx,gm_gx);
            otherwise
                error('Global variable bool_SNTLN_ROOTS must be valid')
        end
        
        [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha_opt,theta_opt);
        
        
    case 'n' % Exclude SNTLN
        
        fx_n = fx_normalised;
        gx_n = gx_normalised;
        
        
        [ux,vx] = GetQuotients(fx_n,gx_n,t,alpha_opt,theta_opt);
        
    otherwise
        error('Global variable')
        
        
end

%%
% Build Sylvester Matrix for normalised, refined coefficients, used in
% comparing singular values.
% EDIT - 23/07/2015 - This has been updated to fx_new rather than fx, since
% this includes SNTLN Struct. pert. Similarly for g.
St_low_rank = BuildSubresultant(fx_n,gx_n,1,alpha_opt,theta_opt);

%% Get the coefficients of the GCD

% Obtain AGCD - Approximate Greatest Common Divisor
switch bool_q
    case 'n'
        % Q is excluded, uw and vw are in form u_{i} \theta \binom
        
        Bi_mk = GetBinomials(m-t)
        Bi_nk = GetBinomials(n-k)
        
        dx = GetGCD(ux,vx,fx_n,gx_n,t,alpha_opt,theta_opt);
    case 'y'
        % Q is included in sylvester matrix,
        switch bool_denom_syl
            case 'y' % denominator included
                dx = GetGCD(ux,vx,fx_n,gx_n,t,alpha_opt,theta_opt);
            case 'n' % denominator excluded
                dx = GetGCD(ux,vx,fx_n,gx_n,t,alpha_opt,theta_opt);
            otherwise
                error('err')
        end
end


switch bool_q
    case 'n'
        % Case n - Q is excluded from sylvester matrix. uv and vw
        % therefore include binomial coefficients which must be removed to
        % obtain coefficients in standard Bernstein Basis.
        
        ux = ux./Bi_mk';
        vx = vx./Bi_nk';
    case 'y'
    otherwise
        error('err')
end

%%
% Apply/Don't Apply structured perturbations to Approximate Polynomial
% Factorisation such that approximation becomes equality.
switch bool_apf
    case 'y'
        % Apply structured perturbations to APF
        switch bool_APF_Roots
            case 'RootSpecificAPF'
                % Use root method which has added constraints.
                
                [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_uk, PostAPF_vk, PostAPF_theta] = ...
                    APF_Roots(fx_n,ux,vx,theta_opt,dx,t);
                
                % Build Post APF_gx
                PostAPF_gx = zeros(m,1);
                for i = 0:1:m-1
                    PostAPF_gx(i+1) = m.*(gm_fx./ gm_gx) .* (PostAPF_fx(i+2) - PostAPF_fx(i+1));
                end
                
                
            case 'StandardAPF'
                [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_uk, PostAPF_vk, PostAPF_theta] = ...
                    APF(fx_n,gx_n,ux,vx,alpha_opt,theta_opt,dx,t);
            otherwise
                error('err')
        end
        
        % update ux,vx,dx values
        dx = PostAPF_dx;
        vx = PostAPF_vk;
        ux = PostAPF_uk;
        
        % Edit 20/07/2015
        fx = PostAPF_fx;
        gx = PostAPF_gx;
        
    case 'n'
        % Dont apply Structured Perturbations to APF
    otherwise
        error('err')
end

%%


f_output = fx;
g_output = gx;
u_output = ux;
v_output = vx;
d_output = dx;




%%
% Assesment of the Sylvester Matrix before processing, post processing, and
% post SNTLN. Before Preprocessing



svd_unproc = svd(St_unproc);
svd_unproc = svd_unproc./norm(svd_unproc);

[~,R] = qr(St_unproc);
R = abs(R);
[R1_rows,~] = size(diag(R));
R1 = R(1:R1_rows,1:R1_rows);
R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);


% After Prerprocessing
svd_preproc = svd(St_preproc);
svd_preproc = svd_preproc./norm(svd_preproc);

[~,R_preproc] = qr(St_preproc);
R_preproc = abs(R_preproc);
[R1_rows,~] = size(diag(R));
% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);
% Get Norms of each row in the matrix R1
R2_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);



% After SNTLN refinement.
svd_low_rank = svd(St_low_rank);
svd_low_rank = svd_low_rank./norm(svd_low_rank);
[~,R_low_rank] = qr(St_low_rank);
R_low_rank = abs(R);
[R1_rows,~] = size(diag(R_low_rank));

% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);

% Get Norms of each row in the matrix R1
R3_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);


% Plot the Singular Values of the Sylvester matrix
switch plot_graphs
    case 1
        figure('name','Singaular values of Sylvester Matrix');
        plot(1:1:length(svd_unproc),log10(svd_unproc),'red-s','DisplayName','Before Preprocessing')
        hold on
        plot(1:1:length(svd_preproc),log10(svd_preproc),'blue-s','DisplayName','After Preprocessing')
        if (bool_sntln == 1)
            plot(1:1:length(svd_low_rank),log10(svd_low_rank),'green-s','DisplayName','With Structured Perturbations')
        end
        legend(gca,'show');
        xlabel('i')
        title('Ordered Singular Values of The Sylvester Matrix S{(f,g)}')
        ylabel('log_{10} Minimal Singular Values ')
        hold off
end







end








