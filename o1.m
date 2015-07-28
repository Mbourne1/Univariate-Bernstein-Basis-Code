function [f_output,g_output,d_output, u_output, v_output,alpha_opt,theta_opt] = ...
    o1(fx,gx)

% This function computes the GCD of two noisy polynomials fx and gx.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                             Inputs:


% ex - (Int) Example Number

% emin - Signal to noise ratio (minimum)

% emax - Signal to noise ratio (maximum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Outputs:

% gcd_output -

% u_output -

% v_output -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                       GLOBAL VARIABLES


% BOOL_APF - (Boolean)
%    1 :- Apply Structured Perturbations to the Approximate Polynomial
%    Factorisation, before obtaining GCD d_{k}.
%    0 :- Don't Include
%    Structured Perturbations.
global bool_apf

% BOOL_Q - (Boolean) Consists of binomial coefficients of coprime
% polynomials.
%   1 :- Include the binomial coefficients from the null space in the
%    Sylvester Matrix.
%   0 :- Exclude the binomial coefficients from the null space in the
%    Sylvester Matrix.
global bool_q

% BOOL_SNTLN - (Boolean)
%    1 :- Include Structured Perturbations in Sylvester Matrix S_{k} before
%    calculating quotients u_{k}, v_{k} and GCD d_{k}.
%    0 :- Don't Include Structured Perturbations.
global bool_sntln

% BOOL_PREPROC - (Boolean)
% It has been shown that the inclusion of preprocessors Geometric mean,
% scaling by alpha, change of independent variable, yield improved results.
%   1 :- Include Preprocessors.
%   0 :- Exclude Preprocessors.
global bool_preproc

% output_format (bool)
% The format of output from file o1.m
%   1 - output u v and d in terms of w (coefficients include theta)
%   0 - output u v and d in terms of x
global bool_denom_syl

% bool_SNTLN_Roots
%   1 - Use SNTLN with g = f' constraint
%   0 - Use SNTLN for two arbitrary polynomials f and g.
global bool_SNTLN_Roots

% bool_APF_Roots
% 1 - Perform APF with g = f' constraint
%   0 - Perform APF for two arbitrary polynomials f and g.
global bool_APF_Roots

% output_format (bool)
% The format of output from file o1.m
%   1 - output u v and d in terms of w (coefficients include theta)
%   0 - output u v and d in terms of x
global output_format;

global bool_plotgraphs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the degree m of polynomial f
m = length(fx) - 1;

% Get the degree n of polynomial g
n = length(gx) - 1;

% GET DEGREE of GCD along with

% subres_unproc - each subresultant S_{k} unprocessed form.

% subres_preproc - each subresultant S_{k} preprocessed,
% includes alpha and theta

% alphas_opt - array of each optimal alpha corresponding to each
% subresultant.

% thetas_opt - array of each optimal theta corresponding to each
% subresultant. gms_fx :- geometric mean of entries of f in each
% subresultant.

% gms_gx :- geometric mean of entries of g in each subresultant.

% Method 1 - 'From Scratch' Build each subresultant from scratch
[t1, subres_unproc1, subres_preproc1, alphas_opt1, thetas_opt1, gms_fx1, gms_gx1] = ...
    GetDegree(fx,gx);

% Method 2 - 'Build Up using S_{k} = A_{k}S_{k-1}B{k}' Build each
% subresultant by matrix manipulation of the previous subresultant.

switch bool_q
    case 0
        % If Q is excluded, we can not build up. since the form of the
        % entries of S_k rely on the inclusion of Q.
        t2 = t1;
    case 1
        [t2, subres_unproc2, subres_preproc2, alphas_opt2, thetas_opt2, gms_fx2, gms_gx2] = ...
            GetDegreeByNewMethod(fx,gx);
end

if t2 == 0
    f_output = fx;
    g_output = gx;
    d_output = 1;
    u_output = 1;
    v_output = 1;
    alpha_opt = 1;
    theta_opt = 1;
    return
end

% Get method. Method 1 - From Scratch Method Method 2 - Build up method.
% Method 2 returns good results for obtaining the degree, but not good
% results for calculating u,v and GCD.


switch bool_q
    case 0
        % if Q is excluded from the sylvester matrix, we were not able to
        % use method 2, so take all values from method 1.
        t = t1;
        subresultant_unproc    = cell2mat(subres_unproc1(t));
        subresultant_preproc   = cell2mat(subres_preproc1(t));
        alpha_opt              = alphas_opt1(t);
        theta_opt              = thetas_opt1(t);
        gm_fx                  = gms_fx1(t);
        gm_gx                  = gms_gx1(t);
        
    case 1
        % if Q was included in the sylvester matrix,
        
        % we chose from method 1 or 2
        method = 1;
        switch method
            case 1
                t = t2;
                subresultant_unproc    = cell2mat(subres_unproc1(t));
                subresultant_preproc   = cell2mat(subres_preproc1(t));
                alpha_opt              = alphas_opt1(t);
                theta_opt              = thetas_opt1(t);
                gm_fx                  = gms_fx1(t);
                gm_gx                  = gms_gx1(t);
            case 2
                t = t2;
                subresultant_unproc    = cell2mat(subres_unproc2(t));
                subresultant_preproc   = cell2mat(subres_preproc2(t));
                alpha_opt              = alphas_opt2(t);
                theta_opt              = thetas_opt2(t);
                gm_fx                  = gms_fx2(t);
                gm_gx                  = gms_gx2(t);
        end
end

format long
fprintf('\n')
fprintf('\nOptimal Value of Alpha: %f \n',alpha_opt);
fprintf('\nOptimal Value of Theta: %f \n',theta_opt);
fprintf('\n')

% Always use the second method ('Build Up') Method to obtain the degree t
% of the GCD.

% If finding the GCD fails, set the degree of the GCD to be 1.
if isempty(t)
    t = 1;
end

% Initialise some useful vectors.
vecm = 0:1:m;
vecn = 0:1:n;
veck = 0:1:t;
vecnk = 0:1:n-t;
vecmk = 0:1:m-t;


% Normalise fx and gx by Geometric mean to obtain fx_n and gx_n Normalise
% by geometric mean obtained from subresultant S_{t}
switch bool_preproc
    case 0 % Exclude preprocessors
        fx_normalised = fx;
        gx_normalised = gx;
        alpha_opt = 1;
        theta_opt = 1;
    case 1 % include preprocessors
        fx_normalised = fx./gm_fx;
        gx_normalised = gx./gm_gx;
end

%
% Build Sylvester Matrix for normalised coefficients, used in comparing
% singular values. Graphed output at end.
Sylvester_1 = Subresultant_BernsteinBasis(fx,gx,1,1,1);


% Obtain f and g in the modified bernstein basis.
fw_n = fx_normalised.*(theta_opt.^(vecm))';
gw_n = gx_normalised.*(theta_opt.^(vecn))';

% Build Sylvester Matrix for normalised coefficients, used in comparing
% singular values. Graphed output at end.
Sylvester_2 = Subresultant_BernsteinBasis(fx_normalised,gx_normalised,theta_opt,alpha_opt,1);

% Get the optimal column of the sylvester matrix to be removed. Where
% removal of the optimal column gives the minmal residual in (Ak x = ck)

[opt_col] = GetOptimalColumn(subresultant_preproc);

% Apply / Don't Apply structured perturbations to Sylvester Matrix NOTE
% SNTLN SHOULD NOT BE USED WITH (METHOD 2 - 'BUILD UP'). SNTLN uses build
% from scratch method.

switch bool_sntln
    case 1
        % Obtain refined values of fx and gx, alpha, theta, and vector x
        % consisting of coefficients of u and v
        
        
        switch bool_SNTLN_Roots
            case 0 % use standard SNTLN
                [fx_new,gx_new,alpha_opt,theta_opt,X] = ...
                    SNTLN(fx_normalised,gx_normalised,alpha_opt,theta_opt,t,opt_col);
                
            case 1 % use roots specific sntln
                [fx_new,gx_new,alpha_opt,theta_opt,X] = ...
                    SNTLN_Roots(fx_normalised,gx_normalised,alpha_opt,theta_opt,t,opt_col,gm_fx,gm_gx);
        end
        % Define X as vecx with the -1 removed in the removed column
        % position.
        format long
        
        vecx =[
            X(1:(opt_col)-1);
            -1;
            X(opt_col:end);
            ];
        
        % Obtain new normalised, polynomials fw and gw in modified
        % Bernstein basis
        fw_n    = fx_new .* (theta_opt.^vecm)';
        gw_n 	= gx_new .* (theta_opt.^vecn)';
        
        vw      = vecx(1:n-t+1);
        uw      = -vecx(n-t+2:end);
        
        
        subresultant_Struct_preprocessed = Subresultant_BernsteinBasis(fx_new,gx_new,theta_opt,alpha_opt,t);
        ck = subresultant_Struct_preprocessed(:,opt_col);
        Ak = subresultant_Struct_preprocessed;
        Ak(:,opt_col) = [];
        
        
        
        %[uw,vw]     = GetQuotients(Ak,ck,n,t,opt_col)
        
        
        
    case 0 % Exclude SNTLN
        
        fx_new = fx_normalised;
        gx_new = gx_normalised;
        
        % Obtain quotient polynomials, by QR Decomposition of the Sylvester
        % Matrix Ax = b, x gives solution vector, containing quotients.
        
        %Edit 11/06/2015
        % suppose we first multiply the coefficient matrix by the thetas contained in uw and vw,
        % such that the vector x_ls only contains u and v in terms of y
        ck = subresultant_preproc(:,opt_col);
        Ak = subresultant_preproc;
        Ak(:,opt_col) = [];
        
        % get \theta^0...theta^{m-t}
        th_mt = zeros(m-t+1,1);
        th_nt = zeros(n-t+1,1);
        
        for i = 0:1:m-t
            th_mt(i+1) = theta_opt^(i);
        end
        
        for i = 0:1:n-t
            th_nt(i+1) = theta_opt^(i);
        end
        
        [uw,vw]     = GetQuotients(Ak,ck,n,t,opt_col);
        ux = uw./th_mt;
        vx = vw./th_nt;

        
        
end


% Build Sylvester Matrix for normalised, refined coefficients, used in
% comparing singular values.
% EDIT - 23/07/2015 - This has been updated to fx_new rather than fx, since
% this includes SNTLN Struct. pert. Similarly for g.
Sylvester_3 = Subresultant_BernsteinBasis(fx_new,gx_new,theta_opt,alpha_opt,1);


% Obtain AGCD - Approximate Greatest Common Divisor
switch bool_q
    case 0 % Q is excluded, uw and vw are in form u_{i} \theta \binom
        Bi_mk = zeros(1,m-t+1);
        for i = 0:1:m-t
            Bi_mk(i+1) = nchoosek(m-t,i);
        end
        Bi_nk = zeros(1,n-t+1);
        for i = 0:1:n-t
            Bi_nk(i+1) = nchoosek(n-t,i);
        end
        dw = GetGCD(uw./Bi_mk',vw./Bi_nk',fw_n,alpha_opt.*gw_n,t);
    case 1 % Q is included in sylvester matrix,
        switch bool_denom_syl
            case 1 % denominator included
                dw = GetGCD(uw,vw,fw_n,alpha_opt.*gw_n,t);
            case 0 % denominator excluded
                dw = GetGCD(uw.*(nchoosek(m+n-t,m-t)),vw.*(nchoosek(m+n-t,n-t)),fw_n,alpha_opt.*gw_n,t);
                %dw = GetGCD(uw,vw,fw_n,alpha_opt.*gw_n,t)
        end
end
% Obtain quotient polynomials and gcd in bernstein basis, with \theta
% removed.

ux = uw./(theta_opt.^vecmk');

vx = vw./(theta_opt.^vecnk');

dx = dw./(theta_opt.^veck');


switch bool_q
    case 0
        % Case 0 - Q is excluded from sylvester matrix. uv and vw
        % therefore include binomial coefficients which must be removed to
        % obtain coefficients in standard Bernstein Basis.
        
        ux = ux./Bi_mk';
        vx = vx./Bi_nk';
end


% Apply/Don't Apply structured perturbations to Approximate Polynomial
% Factorisation such that approximation becomes equality.
switch bool_apf
    case 1
        
        switch bool_APF_Roots
            case 1 % Use root method which has added constraints.
                
                [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_uk, PostAPF_vk, PostAPF_theta] = ...
                    APF_Roots(fx_new,ux,vx,theta_opt,dx,t);
                
                % Build Post APF_gx
                PostAPF_gx = zeros(m,1);
                for i = 0:1:m-1
                   PostAPF_gx(i+1) = m.*(gm_fx./ gm_gx) .* (PostAPF_fx(i+2) - PostAPF_fx(i+1));
                end

                
            case 0
                [PostAPF_fx, PostAPF_gx, PostAPF_dx, PostAPF_uk, PostAPF_vk, PostAPF_theta] = ...
                    APF(fx_new,gx_new,ux,vx,alpha_opt,theta_opt,dx,t);
        end
        
        % update ux,vx,dx values
        dx = PostAPF_dx;
        vx = PostAPF_vk;
        ux = PostAPF_uk;
        
        % Edit 20/07/2015
        fx = PostAPF_fx;
        gx = PostAPF_gx;
        
        % Edit 21/07/2015 - ensure fw and gw are updated after apf.
        fw_n    = PostAPF_fx .* (PostAPF_theta.^vecm)';
        gw_n    = PostAPF_gx .* (PostAPF_theta.^vecn)';
        dw      = PostAPF_dx .* (PostAPF_theta.^veck)';
        uw      = PostAPF_uk .* (PostAPF_theta.^vecmk)';
        vw      = PostAPF_vk .* (PostAPF_theta.^vecnk)';
end



switch output_format
    case 0 % output dx
        f_output = fx;
        g_output = gx;
        u_output = ux;
        v_output = vx;
        d_output = dx;
        
    case 1 % output dw
        f_output = fw_n;
        g_output = gw_n;
        u_output = uw;
        v_output = vw;
        d_output = dw;
        
end




% Assesment of the Sylvester Matrix before processing, post processing, and
% post SNTLN. Before Preprocessing
svd_1 = svd(Sylvester_1);
svd_1 = svd_1./norm(svd_1);

[~,R] = qr(Sylvester_1);
R = abs(R);
[R1_rows,~] = size(diag(R));
% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);
% Get Norms of each row in the matrix R1
R1_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);


% After Prerprocessing
svd_2 = svd(Sylvester_2);
svd_2 = svd_2./norm(svd_2);
[~,R] = qr(Sylvester_2);
R = abs(R);
[R1_rows,~] = size(diag(R));
% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);
% Get Norms of each row in the matrix R1
R2_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);



% After SNTLN refinement.
svd_3 = svd(Sylvester_3);
svd_3 = svd_3./norm(svd_3);
[~,R] = qr(Sylvester_3);
R = abs(R);
[R1_rows,~] = size(diag(R));

% Obtain R1 the top square of the R matrix.
R1 = R(1:R1_rows,1:R1_rows);

% Get Norms of each row in the matrix R1
R3_RowNorm = sqrt(sum(R1.^2,2))./norm(R1);


% Plot the Singular Values of the Sylvester matrix
switch bool_plotgraphs
    case 1
        fig20 = figure(20);
        plot(1:1:length(svd_1),log10(svd_1),'red-s','DisplayName','Before Preprocessing')
        hold on
        plot(1:1:length(svd_2),log10(svd_2),'blue-s','DisplayName','After Preprocessing')
        if (bool_sntln == 1)
            plot(1:1:length(svd_3),log10(svd_3),'green-s','DisplayName','With Structured Perturbations')
        end
        legend(gca,'show');
        xlabel('i')
        title('Ordered Singular Values of The Sylvester Matrix S{(f,g)}')
        ylabel('log_{10} Minimal Singular Values ')
        hold off
end



% Using transposition of the sylvester matrix to obtain coefficients of the
% gcd

GCDMethod = 0;
switch GCDMethod
    case 2 % Use Last non zero row
        a1 = rref(transpose(Sylvester_1));
        a2 = rref(transpose(Sylvester_2));
        a3 = rref(transpose(Sylvester_3));
        a4 = rref(transpose(subresultant_preproc));
        
        
        data = a3;
        data( ~any(data,2), : ) = [];  %rows
        data = data(end,:);
        
        % Get number of cols in the row.
        k = size(data,2);
        
        Bi_k = zeros(1,k);
        for i = 0:1:k-1
            Bi_k(i+1) = nchoosek(k-1,i);
        end
        
        data = data.*Bi_k;
        
        data( :, ~any(data,1) ) = [];  %columns
        data = data./data(1);
        d_output = data';
        
end




end








