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
        
        % Get u(x) and v(x)
        ux_lr = GetWithoutThetas(uw,theta);
        vx_lr = GetWithoutThetas(vw,theta);
        
        alpha_lr = alpha;
        theta_lr = theta;
        
        
        
    case 'Standard SNTLN' % Structured Non-Linear Total Least Norm
        
        % Perform Structured non-linear total least norm SNTLN
        [fx_lr,gx_lr,ux_lr,vx_lr,alpha_lr,theta_lr] = SNTLN(fx,gx,alpha,theta,k,idx_col);
        
        fw = GetWithThetas(fx,theta_lr);
        a_gw = alpha_lr.* GetWithThetas(gx,theta_lr);
        fw_lr = GetWithThetas(fx_lr,theta_lr);
        a_gw_lr = alpha_lr .* GetWithThetas(gx_lr,theta_lr);
        
        S1 = BuildDTQ(fx,gx,k);
        S2 = BuildDTQ(fx_lr,gx_lr,k);
        S3 = BuildDTQ(fw,a_gw,k);
        S4 = BuildDTQ(fw_lr,a_gw_lr,k);
        
        vSingularValues1 = svd(S1);
        vSingularValues2 = svd(S2);
        vSingularValues3 = svd(S3);
        vSingularValues4 = svd(S4);
        
        switch SETTINGS.PLOT_GRAPHS
            case 'y'
                
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
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        alpha_lr = alpha;
        theta_lr = theta;
        
        fw = GetWithThetas(fx,theta);
        gw = GetWithThetas(gx,theta);
        
        % Get quotient polynomials u(x) and v(x)
        [uw,vw] = GetQuotients(fw,alpha.*gw,k);

        % Divide v(w) and u(w) to obtain u(x) and v(x)
        vx_lr = GetWithoutThetas(vw,theta);
        ux_lr = GetWithoutThetas(uw,theta);
        
    otherwise
        error('SETTINGS.LOW_RANK_APPROXIMATION_METHOD must be valid')
end

end