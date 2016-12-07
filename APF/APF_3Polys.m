function [fx_lr, gx_lr, hx_lr, dx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] = ...
    APF_3Polys(fx, gx, hx, ux, vx, wx, alpha, theta, k)
%
% % Inputs
%
% [fx, gx, hx] : Coefficients of polynomial f(x), g(x) and h(x) in the
% Bernstein basis
%
% [ux, vx, wx] : Coefficients of the polynomial u(x), v(x) and w(x) in the 
% Bernstein basis
% 
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
% 
% k : Degree of polynomial d(x)
%
% % Outputs
%
% [fx_lr, gx_lr, hx_lr] : Coefficients of output f(x), g(x) and h(x)
%
% dx_lr : Coefficients of output d(x)
%
% [ux_lr, vx_lr, wx_lr] : Coefficients of output u(x) v(x) and w(x)





global SETTINGS;


switch SETTINGS.APF_METHOD

    case 'Standard APF Nonlinear'
        
        error([mfilename ' : Code Not Yet Developed']);
        
        [fx_lr, gx_lr, hx_lr, dx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] = ...
            APF_Nonlinear_3Polys(fx, gx, hx, ux, vx, wx, alpha, theta, k);
 

    case 'Standard APF Linear'

        error([mfilename ' : Code Not Yet Developed']);
        
        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx,theta);
        a_gw = alpha.*GetWithThetas(gx,theta);
        hw = GetWithThetas(hx,theta);
        
        % Get u(\omega) and v(\omega)
        uw = GetWithThetas(ux,theta);
        vw = GetWithThetas(vx,theta);
        ww = GetWithThetas(wx,theta);
        
        % Get APF and d(\omega)
        [fw_lr, a_gw_lr, dw_lr, uw_lr, vw_lr] = APF_Linear_3Polys(fw, a_gw, hw, uw, vw, ww, k);
        
        % Get f(x) and g(x) after linear APF function
        fx_lr = GetWithoutThetas(fw_lr,theta);
        gx_lr = GetWithoutThetas(a_gw_lr,theta) ./ alpha;
        hx_lr = GetWithoutThetas(hw_lr,theta);
        
        % Get u(x) and v(x) after linear APF function
        ux_lr = GetWithoutThetas(uw_lr,theta);
        vx_lr = GetWithoutThetas(vw_lr,theta);
        wx_lr = GetWithoutThetas(ww_lr,theta);
        
        % Get d(x) after linear APF function
        dx_lr = GetWithoutThetas(dw_lr,theta);
        
        alpha_lr = alpha;
        theta_lr = theta;
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        hx_lr = hx;
        
        ux_lr = ux;
        vx_lr = vx;
        wx_lr = wx;
        
        dx_lr = GetGCDCoefficients_3Polys(ux, vx, wx, fx, gx, hx, k, alpha, theta);
        
        % Get alpha and theta outputs
        alpha_lr = alpha;
        theta_lr = theta;
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename sprintf(' : Error : %s is not a valide APF Method',SETTINGS.APF_METHOD)])
end


end