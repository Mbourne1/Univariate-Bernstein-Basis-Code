function [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = APF_2Polys(fx,gx,ux,vx,alpha,theta,k)
%
% % Inputs
%
% [fx, gx] : Coefficients of polynomials f(x) and g(x) in the Bernstein
% basis
%
% [ux, vx] : Coefficients of the polynomials u(x) and v(x) in the Bernstein basis
% 
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
% 
% k : Degree of polynomial d(x)
%
% % Outputs
%
% [fx_lr, gx_lr] : Coefficients of the polynomials f(x) and g(x)
%
% dx_lr
%
% [ux_lr, vx_lr] : Coefficients of the polynomials u(x) and v(x)





global SETTINGS;


switch SETTINGS.APF_METHOD

    % Standard APF Nonlinear
    % Standard APF Linear
    % None
    
    case 'Standard APF Nonlinear'
        
        [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = ...
            APF_Nonlinear_2Polys(fx,gx,ux,vx,alpha,theta,k);
 

    case 'Standard APF Linear'

        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.*GetWithThetas(gx,theta);
        
        % Get u(\omega) and v(\omega)
        uw = GetWithThetas(ux,theta);
        vw = GetWithThetas(vx,theta);
        
        % Get APF and d(\omega)
        [fw_lr, a_gw_lr, dw_lr, uw_lr, vw_lr] = APF_Linear_2Polys(fw,a_gw,uw,vw,k);
        
        % Get f(x) and g(x) after linear APF function
        fx_lr = GetWithoutThetas(fw_lr,theta);
        gx_lr = GetWithoutThetas(a_gw_lr,theta) ./ alpha;
        
        % Get u(x) and v(x) after linear APF function
        ux_lr = GetWithoutThetas(uw_lr,theta);
        vx_lr = GetWithoutThetas(vw_lr,theta);
        
        % Get d(x) after linear APF function
        dx_lr = GetWithoutThetas(dw_lr,theta);
        
        alpha_lr = alpha;
        theta_lr = theta;
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        ux_lr = ux;
        vx_lr = vx;
        dx_lr = GetGCDCoefficients_2Polys(ux,vx,fx,gx,k,alpha,theta);
        alpha_lr = alpha;
        theta_lr = theta;
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename sprintf(' : Error : %s is not a valid APF Method',SETTINGS.APF_METHOD)])
end


end