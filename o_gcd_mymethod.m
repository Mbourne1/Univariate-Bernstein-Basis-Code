function [fx, gx, dx, ux, vx, alpha, theta, t ] = ...
    o_gcd_mymethod(fx,gx,deg_limits)
% This function computes the GCD d(x) of two noisy polynomials f(x) and g(x).
%
% Inputs:
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


global SETTINGS

% Add relevant Paths


addpath (...
    'APF',...
    'Formatting',...
    'Get GCD Coefficients',...
    'Get GCD Degree',....
    'Get Cofactor Coefficients',...
    'Low Rank Approx');
% Get the degree m of polynomial f
m = GetDegree(fx) ;

% Get degree of GCD by first method

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        figure_name = sprintf('%s : Singular Values',mfilename);
        
        figure('name',figure_name)
        hold on
        title('Singular values of S_{1}')
        plot(log10(svd(BuildSubresultant(fx,gx,1))),'-s');
        hold off
    case 'n'
    otherwise
        error('err');
end

%Get degree by original method - limits
[t, alpha, theta, gm_fx, gm_gx] = Get_GCD_Degree(fx,gx,deg_limits);
LineBreakMedium();

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
[~,idx_col] = GetMinDistance(St_preproc);

% %
% % Get Low rank approximation of the Sylvester matrix S_{t}
% 
[fx_lr,gx_lr,ux_lr,vx_lr,alpha_lr,theta_lr] = LowRankApproximation(fx_n,gx_n,alpha,theta,t,idx_col);


% %
% Build Sylvester Matrix for normalised, refined coefficients, used in
% comparing singular values.
fw = GetWithThetas(fx_lr,theta_lr);
gw = GetWithThetas(gx_lr,theta_lr);
St_low_rank = BuildSubresultant(fw,alpha.*gw,t);

% % 
% Get the coefficients of the GCD by APF or other method.
[fx_alr, gx_alr, dx_alr, ux_alr, vx_alr, alpha_alr, theta_alr] = APF(fx_lr,gx_lr,ux_lr,vx_lr,alpha_lr,theta_lr,t);

fx = fx_alr;
gx = gx_alr;
dx = dx_alr;
ux = ux_alr;
vx = vx_alr;
alpha = alpha_alr;
theta = theta_alr;

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
        plot(log10(vSingVals_St_unproc),'red-s','DisplayName','Before Preprocessing')
        hold on
        plot(log10(vSingVals_St_preproc),'blue-s','DisplayName','After Preprocessing')
        switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
            case {'Standard STLN', 'Standard SNTLN','Root Specific SNTLN'}
                plot(log10(vSingVals_St_LowRank),'green-s','DisplayName','Low Rank Approximation')
            case 'None'
            otherwise
                error('error : Low rank approximation method must be valid')
        end
        legend(gca,'show');
        xlabel('i')
        xlim([1 +inf]);
        title('Ordered Singular Values of The Sylvester Matrix S{(f,g)}')
        ylabel('log_{10} Minimal Singular Values ')
        hold off
    case 'n'
    otherwise
        error('err')
end







end








