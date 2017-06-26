function [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
    LowRankApproximation_2Polys(fx, gx, alpha, theta, k, idx_col)
% Get the low rank approximation of the Sylvester subresultant matrix
% S_{k}(f,g)
%
% Inputs.
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% gx : (Vector) Coefficients of the polynomial g(x)
%
% alpha :  (Float) Optimal value of \alpha
%
% theta : (Float) Optimal value of \theta
%
% k : (Int) Degree of GCD of f(x) and g(x)
%
% idx_col : (Int) Index of optimal column for removal from S_{k}(f,g)
%
% % Outputs
%
% fx_lr : (Vector) Output coefficients of the polynomial f(x)
%
% gx_lr : (Vector) Output coefficients of the polynomial g(x)
%
% ux_lr : (Vector) Output coefficients of the polynomial u(x)
% 
% vx_lr : (Vector) Output coefficients of the polynomial v(x)
%
% alpha_lr : (Float) Optimal value of \alpha
%
% theta_lr : (Float) Optimal value of \theta

global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    
    case 'Standard STLN'
        
        % Get preprocessed polynomials f(\omega) and g(\omega)
        fw = GetWithThetas(fx, theta);
        gw = GetWithThetas(gx, theta);
        
        % Performe STLN to get low rank approximation of S(f,g)
        [fw_lr, a_gw_lr, uw, vw] = STLN(fw, alpha.*gw, k, idx_col);
        
        % Get f(x) and g(x) from low rank approximation.
        fx_lr = GetWithoutThetas(fw_lr, theta);
        gx_lr = GetWithoutThetas(a_gw_lr, theta) ./ alpha;
        
        % Get u(x) and v(x) from low rank approximation
        ux_lr = GetWithoutThetas(uw, theta);
        vx_lr = GetWithoutThetas(vw, theta);
        
        % \alpha and \theta are same as input, as they are unchanged by
        % STLN
        alpha_lr = alpha;
        theta_lr = theta;
        
        %
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, alpha.*gw, fw_lr, a_gw_lr,k);
        
        
        
    case 'Standard SNTLN' % Structured Non-Linear Total Least Norm
        
        % Perform SNTLN to get low rank approximation of S_{k}(f,g)
        [fx_lr, gx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = SNTLN(fx, gx, alpha, theta, k, idx_col);
        
        % Get f(\omega)
        fw = GetWithThetas(fx, theta_lr);
        
        % Get g(\omega)
        a_gw = alpha_lr.* GetWithThetas(gx, theta_lr);
        
        % Get f(\omega) from f(x) from low rank approximation
        fw_lr = GetWithThetas(fx_lr, theta_lr);
        
        % Get g(\omega) from g(x) from low rank approximation
        a_gw_lr = alpha_lr .* GetWithThetas(gx_lr, theta_lr);
        
        % Plot the Singular values of the Sylvester subresultants S_{k}
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, a_gw, fw_lr, a_gw_lr,k)
        
    case 'None'
        
        
        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx, theta);
        gw = GetWithThetas(gx, theta);
        
        % Get quotient polynomials u(\omega) and v(\omega)
        [uw,vw] = GetQuotients_2Polys(fw, alpha.*gw, k);
        
        % f(x) and g(x) are unchanged from input
        fx_lr = fx;
        gx_lr = gx;
        
        % \alpha and \theta are unchanged from input
        alpha_lr = alpha;
        theta_lr = theta;
        
        % Get polynomials u(x) and v(x) from u(\omega) and v(\omega)
        vx_lr = GetWithoutThetas(vw, theta);
        ux_lr = GetWithoutThetas(uw, theta);
        SETTINGS.LOW_RANK_APPROX_REQ_ITE = 0;
    otherwise
        
        error('SETTINGS.LOW_RANK_APPROXIMATION_METHOD must be valid')
        
end

end


function Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr, fw, a_gw, fw_lr, a_gw_lr, k)
%
% % Inputs
%
% fx : (Vector) Coefficients of input polynomial f(x)
%
% gx : (Vector) Coefficients of input polynomial g(x)
%
% fx_lr : (Vector) Coefficients of polynomial f(x) + \delta f(x) from low rank
% approximation method.
%
% gx_lr : (Vector) Coefficients of polynomial g(x) + \delta g(x) from low rank
% approximation method
%
% fw : (Vector) Coefficients of polynomial f(\omega)
%
% gw : (Vector) Coefficients of polynomial g(\omega)
%
% fw_lr : (Vector) 
%
% gw_lr : (Vector)

global SETTINGS

if(SETTINGS.PLOT_GRAPHS)
    
    % Build Sylvester subresultant matrices S_{k} for each of the four
    % pairs of polynomials.
    
    % Unprocessed 
    S1 = BuildDTQ(fx, gx, 1);
    
    % Unprocessed after low rank approximation
    S2 = BuildDTQ(fx_lr, gx_lr, 1);
    
    % Preprocessed
    S3 = BuildDTQ(fw, a_gw, 1);
    
    % Preprocessed after low rank approximation
    S4 = BuildDTQ(fw_lr, a_gw_lr, 1);
    
    % Get singular values for each of the 4 Sylvester subresultant
    % matrices S_{k}
    vSingularValues1 = svd(S1);
    vSingularValues2 = svd(S2);
    vSingularValues3 = svd(S3);
    vSingularValues4 = svd(S4);
    
    % Plot Singular values
    figure_name = sprintf([mfilename ' : Singular Values']);
    figure('name',figure_name)
    plot(log10(vSingularValues1),'-s','DisplayName','fx,gx')
    hold on
    plot(log10(vSingularValues2),'-s','DisplayName','fx_lr,gx_lr')
    plot(log10(vSingularValues3),'-s','DisplayName','fw,a_gw')
    plot(log10(vSingularValues4),'-s','DisplayName','fw_lr,gw_lr')
    legend(gca,'show');
    hold off
    
    
end
end