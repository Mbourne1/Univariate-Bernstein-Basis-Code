function [fx_lr,gx_lr,ux_lr,vx_lr,alpha_lr,theta_lr] = LowRankApproximation(fx,gx,alpha,theta,k,idx_col)
% Get the low rank approximation of the Sylvester subresultant matrix
% S_{k}(f,g)
%
% Inputs.
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x)
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
%
% k : Degree of GCD of f(x) and g(x)
%
% idx_col : Index of optimal column for removal from S_{k}(f,g)
%
% % Outputs
%
% fx_lr :
%
% gx_lr :
%
% ux_lr :
%
% vx_lr :
%
% alpha :
%
% theta :

global SETTINGS

switch SETTINGS.LOW_RANK_APPROXIMATION_METHOD
    case 'Standard STLN'
        
        fw = GetWithThetas(fx,theta);
        gw = GetWithThetas(gx,theta);
        
        % Performe STLN
        [fw_lr,a_gw_lr,uw,vw] = STLN(fw,alpha.*gw,k,idx_col);
        
        % Get f(x) and g(x) from low rank approximation.
        fx_lr = GetWithoutThetas(fw_lr,theta);
        gx_lr = GetWithoutThetas(a_gw_lr,theta) ./ alpha;
        
        % Get u(x) and v(x) from low rank approximation
        ux_lr = GetWithoutThetas(uw,theta);
        vx_lr = GetWithoutThetas(vw,theta);
        
        alpha_lr = alpha;
        theta_lr = theta;
        
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, alpha.*gw, fw_lr, a_gw_lr);
        
        
        
    case 'Standard SNTLN' % Structured Non-Linear Total Least Norm
        
        % Perform Structured non-linear total least norm SNTLN
        [fx_lr,gx_lr,ux_lr,vx_lr,alpha_lr,theta_lr] = SNTLN(fx,gx,alpha,theta,k,idx_col);
        
        % Get f(\omega)
        fw = GetWithThetas(fx,theta_lr);
        
        % Get g(\omega)
        a_gw = alpha_lr.* GetWithThetas(gx,theta_lr);
        
        % Get f(\omega) from f(x) from low rank approximation
        fw_lr = GetWithThetas(fx_lr,theta_lr);
        
        % Get g(\omega) from g(x) from low rank approximation
        a_gw_lr = alpha_lr .* GetWithThetas(gx_lr,theta_lr);
        
        % Plot the Singular values of the Sylvester subresultants S_{k}
        Plot_LowRank_SingularValues(fx, gx, fx_lr, gx_lr,...
            fw, a_gw, fw_lr, a_gw_lr)
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        
        alpha_lr = alpha;
        theta_lr = theta;
        
        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx,theta);
        gw = GetWithThetas(gx,theta);
        
        % Get quotient polynomials u(\omega) and v(\omega)
        [uw,vw] = GetQuotients(fw,alpha.*gw,k);
        
        % Get polynomials u(x) and v(x)
        vx_lr = GetWithoutThetas(vw,theta);
        ux_lr = GetWithoutThetas(uw,theta);
        
    otherwise
        
        error('SETTINGS.LOW_RANK_APPROXIMATION_METHOD must be valid')
        
end

end


function Plot_LowRank_SingularValues(fx,gx,fx_lr,gx_lr,fw,a_gw, fw_lr,a_gw_lr)
%
% % Inputs
%
% fx : Coefficients of input polynomial f(x)
%
% gx : Coefficients of input polynomial g(x)
%
% fx_lr : Coefficients of polynomial f(x) + \delta f(x) from low rank
% approximation method.
%
% gx_lr : Coefficients of polynomial g(x) + \delta g(x) from low rank
% approximation method
%
% fw : Coefficients of polynomial f(\omega)
%
% gw : Coefficients of polynomial g(\omega)
%
% fw_lr :
%
% gw_lr :

global SETTINGS

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        
        % Build Sylvester subresultant matrices S_{k} for each of the four 
        % pairs of polynomials.
        S1 = BuildDTQ(fx,gx,k);
        S2 = BuildDTQ(fx_lr,gx_lr,k);
        S3 = BuildDTQ(fw,a_gw,k);
        S4 = BuildDTQ(fw_lr,a_gw_lr,k);
        
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
        
    case 'n'
    otherwise
        error([mfilename ' : error'])
end
end