function [fx, gx, dx, ux, vx, alpha, theta, t ] = ...
    o_gcd_mymethod(fx,gx,deg_limits)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
%                             Inputs:
%
% fx : Coefficients of the polynomial f(x)
%
% gx : Coefficients of the polynomial g(x)
%
% deg_limits : Upper and lower limits for GCD degree may be defined here
% otherwise set to [0,min(m,n)]
%
%
%
% Outputs:
%
% fx : f(x) + \delta f(x)
%
% gx : g(x) + \delta g(x
%
% dx : The GCD of f(x) + \delta f(x) and g(x) + \delta g(x)
%
% ux : Coefficients of polynomial u(x) where u(x) = f(x)/d(x)
%
% vx : Coefficeints of polynomial v(x) where v(x) = g(x)/d(x)
%
% alpha : Optimal \alpha
%
% theta : Optimal \theta


%
%                       GLOBAL VARIABLES


global SETTINGS

% Get the degree m of polynomial f
m = GetDegree(fx) ;

% Get degree of GCD by first method

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure('name','svd')
        hold on
        plot((svd(BuildSubresultant(fx,gx,1))),'-s');
        hold off
    case 'n'
    otherwise
        error('err');
end

%Get degree by original method - limits
[t2, alpha2, theta2, gm_fx2, gm_gx2] = ...
    GetGCD_Degree(fx,gx,deg_limits);
%display(t2)
LineBreakMedium();
%
% % % Get Degree by new method - no limits
% [t3, alpha3, theta3, gm_fx3, gm_gx3] = ...
%     GetGCD_DegreeByNewMethod(fx,gx);
% display(t3)
% LineBreakMedium();
%
% % Get degree by new method - limits
% [t4, alpha4, theta4, gm_fx4, gm_gx5] = ...
%     GetGCD_DegreeByNewMethod2(fx,gx,deg_limits);
% display(t4)
% LineBreakMedium();

t = t2;
alpha = alpha2;
theta = theta2;
gm_fx = gm_fx2;
gm_gx = gm_gx2;

if t == 0
    % If the two polynomials f(x) and g(x) are coprime, set GCD to be 1,
    fprintf([mfilename ' : ' sprintf('f(x) and g(x) are coprime \n')])
    dx = 1;
    ux = fx;
    vx = gx;
    alpha = 1;
    theta = 1;
    return
end

% If finding the GCD fails, set the degree of the GCD to be 1.
if isempty(t)
    t = 1;
end

% Normalise f(x) and g(x) by Geometric mean to obtain fx_n and gx_n.
% Normalise by geometric mean obtained by entries of f(x) and g(x) in the
% subresultant S_{t}
fx_n = fx./gm_fx;
gx_n = gx./gm_gx;

% Get the Sylvester subresultant of unprocessed f(x) and g(x)
St_unproc = BuildSubresultant(fx,gx,t);

% Get Subresultant of preprocessed f(\theta,\omega) and g(\theta, \omega)
fw = GetWithThetas(fx_n,theta);
gw = GetWithThetas(gx_n,theta);

St_preproc = BuildSubresultant(fw,alpha.*gw,t);

%
% Get the optimal column of the sylvester matrix to be removed. Where
% removal of the optimal column gives the minmal residual in (Ak x = ck)
[opt_col] = GetOptimalColumn(St_preproc);

% %
% % Perform SNTLN
% Apply / Don't Apply structured perturbations.
[fx_n,gx_n,alpha,theta] = LowRankApproximation(fx_n,gx_n,alpha,theta,t,opt_col,gm_fx,gm_gx);


% %
% % Get u(x) and v(x)
% %
% %
fw = GetWithThetas(fx_n,theta);
gw = GetWithThetas(gx_n,theta);

% Get quotient polynomials u(x) and v(x)
[uw,vw] = GetQuotients(fw,alpha.*gw,t);

% Divide v(w) and u(w) to obtain u(x) and v(x)
vx = GetWithoutThetas(vw,theta);
ux = GetWithoutThetas(uw,theta);


% %
% %
% %
% %
% Build Sylvester Matrix for normalised, refined coefficients, used in
% comparing singular values.

fw = GetWithThetas(fx_n,theta);
gw = GetWithThetas(gx_n,theta);

St_low_rank = BuildSubresultant(fw,alpha.*gw,t);

% % Get the coefficients of the GCD
dx = GetGCDCoefficients(ux,vx,fx_n,gx_n,t,alpha,theta);



% %
% Apply/Don't Apply structured perturbations to Approximate Polynomial
% Factorisation such that approximation becomes equality.
switch SETTINGS.APF_METHOD
    case 'Root Specific APF'
        % Use root method which has added constraints.
        
        [APF_fx, APF_gx, APF_dx, APF_uk, APF_vk, APF_theta] = ...
            APF_Roots(fx_n,ux,vx,theta,dx,t);
        
        % Build Post APF_gx
        APF_gx = zeros(m,1);
        
        for i = 0:1:m-1
            APF_gx(i+1) = m.*(gm_fx./ gm_gx) .* (APF_fx(i+2) - APF_fx(i+1));
        end
        
        % update ux,vx,dx values
        dx = APF_dx;
        vx = APF_vk;
        ux = APF_uk;
        
        % Edit 20/07/2015
        fx = APF_fx;
        gx = APF_gx;
        
    case 'Standard APF'
        [APF_fx, APF_gx, APF_dx, APF_uk, APF_vk, APF_theta] = ...
            APF(fx_n,gx_n,ux,vx,alpha,theta,dx,t);
        
        % update ux,vx,dx values
        dx = APF_dx;
        vx = APF_vk;
        ux = APF_uk;
        
        % Edit 20/07/2015
        fx = APF_fx;
        gx = APF_gx;
    case 'None'
        % Do Nothing
    otherwise
        error('SETTINGS.APF_METHOD IS either Standard APF, Root Specific APF or None')
end






% %
% Assesment of the Sylvester Matrix before processing, post processing, and
% post SNTLN. Before Preprocessing


% Get singular values of Sylvester matrix of noisy input polynomials
vSingVals_St_unproc = svd(St_unproc);
vSingVals_St_unproc = vSingVals_St_unproc./norm(vSingVals_St_unproc);

% Get Singular values of Sylvester matrix of preprocessed polynomials
vSingVals_St_preproc = svd(St_preproc);
vSingVals_St_preproc = vSingVals_St_preproc./norm(vSingVals_St_preproc);

% Get Singular values of Sylvester matrix of f+ \delta f and g+\delta g
vSingVals_St_LowRank = svd(St_low_rank);
vSingVals_St_LowRank = vSingVals_St_LowRank./norm(vSingVals_St_LowRank);


% Plot the Singular Values of the Sylvester matrix
switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure('name','Singaular values of Sylvester Matrix');
        plot(1:1:length(vSingVals_St_unproc),log10(vSingVals_St_unproc),'red-s','DisplayName','Before Preprocessing')
        hold on
        plot(1:1:length(vSingVals_St_preproc),log10(vSingVals_St_preproc),'blue-s','DisplayName','After Preprocessing')
        switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
            case {'Standard STLN', 'Standard SNTLN','Root Specific SNTLN'}
                plot(1:1:length(vSingVals_St_LowRank),log10(vSingVals_St_LowRank),'green-s','DisplayName','Low Rank Approximation')
            case 'None'
            otherwise
                error('error : Low rank approximation method must be valid')
        end
        legend(gca,'show');
        xlabel('i')
        title('Ordered Singular Values of The Sylvester Matrix S{(f,g)}')
        ylabel('log_{10} Minimal Singular Values ')
        hold off
    case 'n'
    otherwise
        error('err')
end







end








