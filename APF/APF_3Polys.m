function [fx_lr, gx_lr, hx_lr, dx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, beta_lr, theta_lr] = ...
    APF_3Polys(fx, gx, hx, ux, vx, wx, alpha, beta, theta, k)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% gx : (Vector) Coefficients of polynomial g(x)
%
% hx : (Vector) Coefficients of polynomial h(x)
%
% ux : (Vector) Coefficients of polynomial u(x)
%
% vx : (Vector) Coefficients of polynomial v(x)
%
% wx : (Vector) Coefficients of the polynomial w(x)
%
% alpha : (Float) Optimal value of \alpha
%
% theta : (Float) Optimal value of \theta
% 
% k : (Int) Degree of polynomial d(x)
%
% % Outputs
%
% fx_lr : (Vector) Coefficients of output f(x)
%
% gx_lr : (Vector) Coefficients of output g(x)
%
% hx_lr : (Vector) Coefficients of output h(x)
%
% dx_lr : (Vector) Coefficients of output d(x)
%
% ux_lr : (Vector) Coefficients of output u(x)
%
% vx_lr : (Vector) Coefficients of output v(x)
%
% wx_lr : (Vector) Coefficients of output w(x)





global SETTINGS;


switch SETTINGS.APF_METHOD

    case 'Standard APF Nonlinear'
        
        error([mfilename ' : Code Not Yet Developed']);
        
        % This function does not exist
        
        [fx_lr, gx_lr, hx_lr, dx_lr, ux_lr, vx_lr, wx_lr, alpha_lr, theta_lr] = ...
            APF_Nonlinear_3Polys(fx, gx, hx, ux, vx, wx, alpha, beta, theta, k);
 

    case 'Standard APF Linear'

        error([mfilename ' : Code Not Yet Developed']);
        
        % Get f(\omega) and g(\omega)
        fw = GetWithThetas(fx, theta);
        a_gw = alpha.*GetWithThetas(gx, theta);
        b_hw = GetWithThetas(hx, theta);
        
        % Get u(\omega) and v(\omega)
        uw = GetWithThetas(ux, theta);
        vw = GetWithThetas(vx, theta);
        ww = GetWithThetas(wx, theta);
        
        % Get APF and d(\omega)
        [fw_lr, a_gw_lr, dw_lr, uw_lr, vw_lr] = APF_Linear_3Polys(fw, a_gw, b_hw, uw, vw, ww, k);
        
        % Get f(x) and g(x) after linear APF function
        fx_lr = GetWithoutThetas(fw_lr, theta);
        gx_lr = GetWithoutThetas(a_gw_lr, theta) ./ alpha;
        hx_lr = GetWithoutThetas(hw_lr, theta) ./ beta;
        
        % Get u(x) and v(x) after linear APF function
        ux_lr = GetWithoutThetas(uw_lr, theta);
        vx_lr = GetWithoutThetas(vw_lr, theta);
        wx_lr = GetWithoutThetas(ww_lr, theta);
        
        % Get d(x) after linear APF function
        dx_lr = GetWithoutThetas(dw_lr, theta);
        
        alpha_lr = alpha;
        beta_lr = beta;
        theta_lr = theta;
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        hx_lr = hx;
        
        ux_lr = ux;
        vx_lr = vx;
        wx_lr = wx;
        
        dx_lr = GetGCDCoefficients_3Polys(ux, vx, wx, fx, gx, hx, k, alpha, beta, theta);
        
        % Get alpha and theta outputs
        alpha_lr = alpha;
        beta_lr = beta;
        theta_lr = theta;
        
        SETTINGS.APF_REQ_ITE = 0;
        
    otherwise
        error([mfilename sprintf(' : Error : %s is not a valide APF Method',SETTINGS.APF_METHOD)])
end


end