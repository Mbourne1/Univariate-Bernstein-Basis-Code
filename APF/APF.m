function [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = APF(fx,gx,ux,vx,alpha,theta,k)
%
% % Inputs
%
% fx : Coefficients of polynomial f(x)
%
% gx : Coefficients of polynomial g(x) in the Bernstein basis
%
% ux : Coefficients of the polynomial u(x) in the Bernstein basis
% 
% vx : Coefficients of the polynomial v(x) in the Bernstein basis
%
% alpha : Optimal value of \alpha
%
% theta : Optimal value of \theta
% 
% k : Degree of polynomial d(x)
%
% % Outputs
%
% fx_lr
%
% gx_lr 
%
% dx_lr
%
% ux_lr
%
% vx_lr





global SETTINGS;


switch SETTINGS.APF_METHOD

    case 'Standard APF NonLinear'
        
        [fx_lr, gx_lr, dx_lr, ux_lr, vx_lr, alpha_lr, theta_lr] = APF_NonLinear(fx,gx,ux,vx,alpha,theta,k);
        
        
        
        
    case 'Standard APF Linear'
        
        error('Not Developed')
        APF_Linear()
        
    case 'None'
        
        fx_lr = fx;
        gx_lr = gx;
        ux_lr = ux;
        vx_lr = vx;
        dx_lr = GetGCDCoefficients(ux,vx,fx,gx,k,alpha,theta);
        alpha_lr = alpha;
        theta_lr = theta;
    
end


end